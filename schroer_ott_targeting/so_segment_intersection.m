function hit = so_segment_intersection(a, b, c, d, tol)
%SO_SEGMENT_INTERSECTION Exact 2-D segment intersection parameters.
if nargin < 5
    tol = 1e-12;
end
r = b - a;
s = d - c;
den = cross2(r, s);
qmp = c - a;
hit.success = false;
hit.alpha = NaN;
hit.beta = NaN;
hit.point = [NaN; NaN];

if abs(den) <= tol
    if abs(cross2(qmp, r)) > tol
        return;
    end
    rr = dot(r, r);
    if rr <= tol
        return;
    end
    t0 = dot(c - a, r) / rr;
    t1 = dot(d - a, r) / rr;
    lo = max(0, min(t0, t1));
    hi = min(1, max(t0, t1));
    if lo <= hi + tol
        hit.success = true;
        hit.alpha = max(0, min(1, 0.5 * (lo + hi)));
        point = a + hit.alpha * r;
        ss = dot(s, s);
        if ss > tol
            hit.beta = max(0, min(1, dot(point - c, s) / ss));
        else
            hit.beta = 0;
        end
        hit.point = point;
    end
    return;
end

alpha = cross2(qmp, s) / den;
beta = cross2(qmp, r) / den;
if alpha >= -tol && alpha <= 1 + tol && beta >= -tol && beta <= 1 + tol
    hit.success = true;
    hit.alpha = max(0, min(1, alpha));
    hit.beta = max(0, min(1, beta));
    hit.point = a + hit.alpha * r;
end
end

function z = cross2(u, v)
z = u(1) * v(2) - u(2) * v(1);
end

