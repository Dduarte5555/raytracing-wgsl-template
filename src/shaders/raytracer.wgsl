const THREAD_COUNT = 16;
const RAY_TMIN = 0.0001;
const RAY_TMAX = 100.0;
const PI = 3.1415927f;
const FRAC_1_PI = 0.31830987f;
const FRAC_2_PI = 1.5707964f;

@group(0) @binding(0)  
  var<storage, read_write> fb : array<vec4f>;

@group(0) @binding(1)
  var<storage, read_write> rtfb : array<vec4f>;

@group(1) @binding(0)
  var<storage, read_write> uniforms : array<f32>;

@group(2) @binding(0)
  var<storage, read_write> spheresb : array<sphere>;

@group(2) @binding(1)
  var<storage, read_write> quadsb : array<quad>;

@group(2) @binding(2)
  var<storage, read_write> boxesb : array<box>;

@group(2) @binding(3)
  var<storage, read_write> trianglesb : array<triangle>;

@group(2) @binding(4)
  var<storage, read_write> meshb : array<mesh>;

struct ray {
  origin : vec3f,
  direction : vec3f,
};

struct sphere {
  transform : vec4f,
  color : vec4f,
  material : vec4f,
};

struct quad {
  Q : vec4f,
  u : vec4f,
  v : vec4f,
  color : vec4f,
  material : vec4f,
};

struct box {
  center : vec4f,
  radius : vec4f,
  rotation: vec4f,
  color : vec4f,
  material : vec4f,
};

struct triangle {
  v0 : vec4f,
  v1 : vec4f,
  v2 : vec4f,
};

struct mesh {
  transform : vec4f,
  scale : vec4f,
  rotation : vec4f,
  color : vec4f,
  material : vec4f,
  min : vec4f,
  max : vec4f,
  show_bb : f32,
  start : f32,
  end : f32,
};

struct material_behaviour {
  scatter : bool,
  direction : vec3f,
};

struct camera {
  origin : vec3f,
  lower_left_corner : vec3f,
  horizontal : vec3f,
  vertical : vec3f,
  u : vec3f,
  v : vec3f,
  w : vec3f,
  lens_radius : f32,
};

struct hit_record {
  t : f32,
  p : vec3f,
  normal : vec3f,
  object_color : vec4f,
  object_material : vec4f,
  frontface : bool,
  hit_anything : bool,
};

fn ray_at(r: ray, t: f32) -> vec3f
{
  return r.origin + t * r.direction;
}

fn get_ray(cam: camera, uv: vec2f, rng_state: ptr<function, u32>) -> ray
{
  var rd = cam.lens_radius * rng_next_vec3_in_unit_disk(rng_state);
  var offset = cam.u * rd.x + cam.v * rd.y;
  return ray(cam.origin + offset, normalize(cam.lower_left_corner + uv.x * cam.horizontal + uv.y * cam.vertical - cam.origin - offset));
}

fn get_camera(lookfrom: vec3f, lookat: vec3f, vup: vec3f, vfov: f32, aspect_ratio: f32, aperture: f32, focus_dist: f32) -> camera
{
  var camera = camera();
  camera.lens_radius = aperture / 2.0;

  var theta = degrees_to_radians(vfov);
  var h = tan(theta / 2.0);
  var w = aspect_ratio * h;

  camera.origin = lookfrom;
  camera.w = normalize(lookfrom - lookat);
  camera.u = normalize(cross(vup, camera.w));
  camera.v = cross(camera.u, camera.w);

  camera.lower_left_corner = camera.origin - w * focus_dist * camera.u - h * focus_dist * camera.v - focus_dist * camera.w;
  camera.horizontal = 2.0 * w * focus_dist * camera.u;
  camera.vertical = 2.0 * h * focus_dist * camera.v;

  return camera;
}


fn envoriment_color(direction: vec3f, color1: vec3f, color2: vec3f) -> vec3f
{
  var unit_direction = normalize(direction);
  var t = 0.5 * (unit_direction.y + 1.0);
  var col = (1.0 - t) * color1 + t * color2;

  var sun_direction = normalize(vec3(uniforms[13], uniforms[14], uniforms[15]));
  var sun_color = int_to_rgb(i32(uniforms[17]));
  var sun_intensity = uniforms[16];
  var sun_size = uniforms[18];

  var sun = clamp(dot(sun_direction, unit_direction), 0.0, 1.0);
  col += sun_color * max(0, (pow(sun, sun_size) * sun_intensity));

  return col;
}

fn check_ray_collision(r: ray, max: f32) -> hit_record {
    var spheresCount = i32(uniforms[19]);
    var quadsCount = i32(uniforms[20]);
    var boxesCount = i32(uniforms[21]);
    var trianglesCount = i32(uniforms[22]);
    var meshCount = i32(uniforms[27]);

    var record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
    var closest_t = max;

    for (var i = 0; i < spheresCount; i = i + 1) {
        var sphere = spheresb[i];
        var temp_record = hit_record();

        if (hit_sphere(sphere.transform.xyz, sphere.transform.w, r, &temp_record, closest_t)) {
            closest_t = temp_record.t;
            record = temp_record;
            record.object_color = sphere.color;
            record.object_material = sphere.material;
        }
    }

    for (var i = 0; i < trianglesCount; i = i + 1) {
        var triangle = trianglesb[i];
        var temp_record = hit_record();

        if (hit_triangle(r, triangle.v0.xyz, triangle.v1.xyz, triangle.v2.xyz, &temp_record, closest_t)) {
            closest_t = temp_record.t;
            record = temp_record;
        }
    }

    for (var i = 0; i < boxesCount; i = i + 1) {
        var box = boxesb[i];
        if (AABB_intersect(r, box.center.xyz - box.radius.xyz, box.center.xyz + box.radius.xyz)) {
            record.hit_anything = true;
        }
    }

    if (closest_t < max) {
        record.hit_anything = true;
    }

    return record;
}

fn lambertian(normal: vec3f, absorption: f32, random_sphere: vec3f, rng_state: ptr<function, u32>) -> material_behaviour {
    var scatter_direction = normal + rng_next_vec3_in_unit_sphere(rng_state);

    // Evitar a direção zero (caso o vetor aleatório cancele a normal)
    if (all(scatter_direction == vec3f(0.0))) {
        scatter_direction = normal;
    }

    return material_behaviour(true, normalize(scatter_direction));
}

fn metal(normal: vec3f, direction: vec3f, fuzz: f32, random_sphere: vec3f) -> material_behaviour {
    var reflected = reflect(normalize(direction), normalize(normal));
    var scattered = reflected + fuzz * random_sphere;

    if (dot(scattered, normal) <= 0.0) {
        return material_behaviour(false, vec3f(0.0));
    }

    return material_behaviour(true, normalize(scattered));
}

fn refract(uv: vec3f, n: vec3f, eta: f32) -> vec3f {
    var cos_theta = dot(-uv, n);
    var r_out_perp = eta * (uv + cos_theta * n);
    var r_out_parallel = -sqrt(abs(1.0 - dot(r_out_perp, r_out_perp))) * n;
    return r_out_perp + r_out_parallel;
}

fn reflectance(cosine: f32, ref_idx: f32) -> f32 {
    var r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) * pow(1.0 - cosine, 5.0);
}

fn dielectric(normal : vec3f, r_direction: vec3f, refraction_index: f32, frontface: bool, random_sphere: vec3f, fuzz: f32, rng_state: ptr<function, u32>) -> material_behaviour
{  
    var refract_ratio = select(refraction_index, 1.0 / refraction_index, frontface);
    var unit_direction = normalize(r_direction);

    var cos_theta = min(dot(-unit_direction, normal), 1.0);
    var sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    var cannot_refract = refract_ratio * sin_theta > 1.0;

    var direction: vec3f;
    if (cannot_refract || reflectance(cos_theta, refract_ratio) > rng_next_float(rng_state)) {
        direction = reflect(unit_direction, normal);
    } else {
        direction = refract(unit_direction, normal, refract_ratio);
    }

    direction += fuzz * random_sphere;

    return material_behaviour(true, normalize(direction));
}

fn emissive(color: vec3f, intensity: f32) -> material_behaviour {
    // Garantir que a cor e a intensidade estão dentro de limites aceitáveis
    var clamped_intensity = clamp(intensity, 0.0, 10.0); // Limite máximo opcional
    var clamped_color = clamp(color, vec3f(0.0), vec3f(1.0)); // Garante que as cores estão entre 0 e 1

    // Retornar o comportamento com a luz emitida
    return material_behaviour(false, clamped_color * clamped_intensity);
}

fn trace(r: ray, rng_state: ptr<function, u32>) -> vec3f {
    var maxbounces = i32(uniforms[2]);
    var light = vec3f(0.0);
    var color = vec3f(1.0);
    var r_ = r;

    var backgroundcolor1 = int_to_rgb(i32(uniforms[11]));
    var backgroundcolor2 = int_to_rgb(i32(uniforms[12]));
    var behaviour = material_behaviour(true, vec3f(0.0));

    for (var j = 0; j < maxbounces; j = j + 1) {
        var record = check_ray_collision(r_, RAY_TMAX);

        // Se não houver interseção, retorna cor do fundo
        if (!record.hit_anything) {
            light += color * envoriment_color(r_.direction, backgroundcolor1, backgroundcolor2);
            break;
        }

        // Processar comportamento do material no ponto de interseção
        if (record.object_material.x == 0.0) {
            behaviour = lambertian(record.normal, 1.0, vec3f(0.0), rng_state);
        } else if (record.object_material.x == 1.0) {
            var random_sphere = rng_next_vec3_in_unit_sphere(rng_state);
            behaviour = metal(record.normal, r_.direction, record.object_material.y, random_sphere);
            color *= vec3f(1.0); // Não aplica a cor diretamente para metais
        } else if (record.object_material.x == 2.0) {
            behaviour = dielectric(
                record.normal,
                r_.direction,
                record.object_material.y,
                record.frontface,
                vec3f(0.0),
                record.object_material.z,
                rng_state
            );
        } else if (record.object_material.x == 3.0) {
             behaviour = emissive(record.object_color.rgb, record.object_material.y);
    
              // Limitar intensidade da emissão (opcional)
              var intensity = min(record.object_material.y, 10.0);
              
              // Adicionar luz emitida diretamente
              light += intensity * behaviour.direction;

              // Interromper o traçamento, já que não há espalhamento
              break;
        }

        // Atualizar o raio para o próximo bounce
        r_ = ray(record.p, behaviour.direction);
        color *= record.object_color.rgb; // Multiplica para difusos e dielétricos

        // Se o material não espalha (absorve), interrompe o loop
        if (!behaviour.scatter) {
            break;
        }
    }

    return light;
}

fn accumulate_framebuffer(index: u32, color_out: vec4f) -> vec4f {
    return mix(rtfb[index], color_out, 0.5);
}

@compute @workgroup_size(THREAD_COUNT, THREAD_COUNT, 1)
fn render(@builtin(global_invocation_id) id: vec3u) {
    var rez = u32(uniforms[1]);
    var time = u32(uniforms[0]);
    var rng_state = init_rng(vec2(id.x, id.y), vec2(rez), time);
    var fragCoord = vec2f(f32(id.x), f32(id.y));
    var uv = (fragCoord + sample_square(&rng_state)) / vec2(f32(rez));

    var lookfrom = vec3(uniforms[7], uniforms[8], uniforms[9]);
    var lookat = vec3(uniforms[23], uniforms[24], uniforms[25]);
    var cam = get_camera(
        lookfrom,
        lookat,
        vec3(0.0, 1.0, 0.0),
        uniforms[10],
        1.0,
        uniforms[6],
        uniforms[5]
    );

    var samples_per_pixel = i32(uniforms[4]);
    var color = vec3f(0.0);

    for (var s = 0; s < samples_per_pixel; s = s + 1) {
        var ray = get_ray(cam, uv, &rng_state);
        var sample_color = trace(ray, &rng_state);
        color = color + sample_color;
    }

    color = color / f32(samples_per_pixel);
    var color_out = vec4(linear_to_gamma(color), 1.0);
    var map_fb = mapfb(id.xy, f32(rez)); // mapfb ainda retorna i32

    var should_accumulate = uniforms[3];
    if (should_accumulate > 0.0) {
        color_out = accumulate_framebuffer(u32(map_fb), color_out); // Conversão para u32
    }

    rtfb[u32(map_fb)] = color_out;
    fb[u32(map_fb)] = color_out;
}