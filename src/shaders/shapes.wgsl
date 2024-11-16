fn hit_sphere(center: vec3f, radius: f32, r: ray, record: ptr<function, hit_record>, max: f32) -> bool {
    var oc = r.origin - center;
    var a = dot(r.direction, r.direction);
    var b = 2.0 * dot(oc, r.direction);
    var c = dot(oc, oc) - radius * radius;
    var discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0.0) {
        return false;
    }

    var sqrtd = sqrt(discriminant);

    // Verificar as possíveis raízes
    var t1 = (-b - sqrtd) / (2.0 * a);
    if (t1 > RAY_TMIN && t1 < max) {
        (*record).t = t1;
        (*record).p = ray_at(r, t1);
        (*record).normal = normalize((*record).p - center);
        return true;
    }

    var t2 = (-b + sqrtd) / (2.0 * a);
    if (t2 > RAY_TMIN && t2 < max) {
        (*record).t = t2;
        (*record).p = ray_at(r, t2);
        (*record).normal = normalize((*record).p - center);
        return true;
    }

    return false;
}

fn hit_quad(r: ray, Q: vec4f, u: vec4f, v: vec4f, record: ptr<function, hit_record>, max: f32)
{
  var n = cross(u.xyz, v.xyz);
  var normal = normalize(n);
  var D = dot(normal, Q.xyz);
  var w = n / dot(n.xyz, n.xyz);

  var denom = dot(normal, r.direction);
  if (abs(denom) < 0.0001)
  {
    record.hit_anything = false;
    return;
  }

  var t = (D - dot(normal, r.origin)) / denom;
  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  var intersection = ray_at(r, t);
  var planar_hitpt_vector = intersection - Q.xyz;
  var alpha = dot(w, cross(planar_hitpt_vector, v.xyz));
  var beta = dot(w, cross(u.xyz, planar_hitpt_vector));

  if (alpha < 0.0 || alpha > 1.0 || beta < 0.0 || beta > 1.0)
  {
    record.hit_anything = false;
    return;
  }

  if (dot(normal, r.direction) > 0.0)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = intersection;
  record.normal = normal;
  record.hit_anything = true;
}

fn hit_triangle(r: ray, v0: vec3f, v1: vec3f, v2: vec3f, record: ptr<function, hit_record>, max: f32) -> bool {
    var v1v0 = v1 - v0;
    var v2v0 = v2 - v0;
    var rov0 = r.origin - v0;

    var n = cross(v1v0, v2v0);
    var q = cross(rov0, r.direction);

    var d = dot(r.direction, n);
    if (abs(d) < 1e-8) {
        return false; // Evita divisão por zero (raio paralelo ao triângulo)
    }

    var inv_d = 1.0 / d;

    var u = inv_d * dot(-q, v2v0);
    if (u < 0.0 || u > 1.0) {
        return false; // Fora do triângulo
    }

    var v = inv_d * dot(q, v1v0);
    if (v < 0.0 || (u + v) > 1.0) {
        return false; // Fora do triângulo
    }

    var t = inv_d * dot(-n, rov0);
    if (t < RAY_TMIN || t > max) {
        return false; // Fora do intervalo válido
    }

    (*record).t = t;
    (*record).p = ray_at(r, t);
    (*record).normal = normalize(n);
    (*record).hit_anything = true;

    return true;
}

fn hit_box(r: ray, center: vec3f, rad: vec3f, record: ptr<function, hit_record>, t_max: f32)
{
  var m = 1.0 / r.direction;
  var n = m * (r.origin - center);
  var k = abs(m) * rad;

  var t1 = -n - k;
  var t2 = -n + k;

  var tN = max(max(t1.x, t1.y), t1.z);
  var tF = min(min(t2.x, t2.y), t2.z);

  if (tN > tF || tF < 0.0)
  {
    record.hit_anything = false;
    return;
  }

  var t = tN;
  if (t < RAY_TMIN || t > t_max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = ray_at(r, t);
  record.normal = -sign(r.direction) * step(t1.yzx, t1.xyz) * step(t1.zxy, t1.xyz);
  record.hit_anything = true;

  return;
}