typedef struct IntersectionProblemParameters {
    double ax;
    double ay;
    double bx;
    double by;
    double mx;
    double my;
    double r;
} IntersectionProblemParameters_t;


double intersection_alpha(IntersectionProblemParameters_t p)
{
    return sqr(p.by - p.ay) + sqr(p.bx - p.ax);
}

double intersection_beta(IntersectionProblemParameters_t p)
{
    return 2 * p.ay * (p.by - p.ay) 
        + 2 * p.ax * (p.bx - p.ax) 
        - 2 * p.my * (p.by - p.ay) 
        - 2 * p.mx * (p.bx - p.ax);
}

double intersection_gamma(IntersectionProblemParameters_t p)
{
    return p.ay * p.ay - 2 * p.ay * p.my + p.my * p.my 
        - p.r * p.r + p.ax * p.ax - 2 * p.ax * p.mx + p.mx * p.mx;
}

double intersection_discriminant(IntersectionProblemParameters_t p)
{
    return sqr(intersection_beta(p)) - 4 * intersection_alpha(p) * intersection_gamma(p);
}

double intersection_s(IntersectionProblemParameters_t p, int sign)
{
    double scale_parameter = (-intersection_beta(p) + 
        sign * my_sqrt(intersection_discriminant(p))) / 2 / intersection_alpha(p);
    if (scale_parameter < 0) scale_parameter = my_nan();
    return scale_parameter;
}

double intersection_s1(IntersectionProblemParameters_t p)
{
    return intersection_s(p, -1);
}

double intersection_s2(IntersectionProblemParameters_t p)
{
    return intersection_s(p, +1);
}

double intersection_x1(IntersectionProblemParameters_t p)
{
    return p.ax + (p.bx - p.ax) * intersection_s1(p);
}

double intersection_x2(IntersectionProblemParameters_t p)
{
    return p.ax + (p.bx - p.ax) * intersection_s2(p);
}

double intersection_y1(IntersectionProblemParameters_t p)
{
    return p.ay + (p.by - p.ay) * intersection_s1(p);
}

double intersection_y2(IntersectionProblemParameters_t p)
{
    return p.ay + (p.by - p.ay) * intersection_s2(p);
}

int number_of_intersections(IntersectionProblemParameters_t p)
{
    double d = intersection_discriminant(p);
    if (d < 0) return 0;
    else if (d == 0) return 1;
    else if (d > 0) {
        // Both intersection points behind the trajectory starting point A:
        if (my_is_nan(intersection_s2(p))) return 0; 
        // One intersection point behind the trajectory starting point A:
        if (my_is_nan(intersection_s1(p))) return 1;
    }
    // Both intersection points on the positive trajectory:
    return 2;
}

double squared_distance_from_center(double X, double Y, double MX, double MY)
{
    return (sqr(MX - X) + sqr(MY - Y));
}

bool intersecting_trajectory_starts_inside(IntersectionProblemParameters_t p)
{
    return (squared_distance_from_center(p.ax, p.ay, p.mx, p.my) < sqr(p.r));
}
