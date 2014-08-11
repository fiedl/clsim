typedef struct IntersectionProblemParameters {
    floating_t ax;
    floating_t ay;
    floating_t bx;
    floating_t by;
    floating_t mx;
    floating_t my;
    floating_t r;
} IntersectionProblemParameters_t;


floating_t intersection_alpha(IntersectionProblemParameters_t p)
{
    return sqr(p.by - p.ay) + sqr(p.bx - p.ax);
}

floating_t intersection_beta(IntersectionProblemParameters_t p)
{
    return 2 * p.ay * (p.by - p.ay) 
        + 2 * p.ax * (p.bx - p.ax) 
        - 2 * p.my * (p.by - p.ay) 
        - 2 * p.mx * (p.bx - p.ax);
}

floating_t intersection_gamma(IntersectionProblemParameters_t p)
{
    return p.ay * p.ay - 2 * p.ay * p.my + p.my * p.my 
        - p.r * p.r + p.ax * p.ax - 2 * p.ax * p.mx + p.mx * p.mx;
}

floating_t intersection_discriminant(IntersectionProblemParameters_t p)
{
    return sqr(intersection_beta(p)) - 4 * intersection_alpha(p) * intersection_gamma(p);
}

floating_t intersection_s(IntersectionProblemParameters_t p, int sign)
{
    floating_t scale_parameter = (-intersection_beta(p) + 
        sign * my_sqrt(intersection_discriminant(p))) / 2 / intersection_alpha(p);
    if (( scale_parameter < 0.0 ) || ( scale_parameter > 1.0 )) scale_parameter = my_nan();
    return scale_parameter;
}

floating_t intersection_s1(IntersectionProblemParameters_t p)
{
    return intersection_s(p, -1);
}

floating_t intersection_s2(IntersectionProblemParameters_t p)
{
    return intersection_s(p, +1);
}

floating_t intersection_x1(IntersectionProblemParameters_t p)
{
    return p.ax + (p.bx - p.ax) * intersection_s1(p);
}

floating_t intersection_x2(IntersectionProblemParameters_t p)
{
    return p.ax + (p.bx - p.ax) * intersection_s2(p);
}

floating_t intersection_y1(IntersectionProblemParameters_t p)
{
    return p.ay + (p.by - p.ay) * intersection_s1(p);
}

floating_t intersection_y2(IntersectionProblemParameters_t p)
{
    return p.ay + (p.by - p.ay) * intersection_s2(p);
}

int number_of_intersections(IntersectionProblemParameters_t p)
{
    floating_t d = intersection_discriminant(p);

    // if (my_is_nan(d))
    // {
    //     printf("DISCRIMINANT ERROR: d = %f\n", d);
    //     printf(" -> alpha = %f\n", intersection_alpha(p));
    //     printf(" -> beta  = %f\n", intersection_beta(p));
    //     printf(" -> gamma = %f\n", intersection_gamma(p));
    //     printf(" -> start at: (%f, %f)\n", p.ax, p.ay);
    //     printf(" -> end at:   (%f, %f)\n", p.bx, p.by);
    //     printf(" -> circle:   (%f, %f), %f\n", p.mx, p.my, p.r);
    // }

    if (d < 0) return 0;
    else if (d == 0) return 1;
    else if (d > 0) {
        // Both intersection points behind the trajectory starting point A:
        if (my_is_nan(intersection_s2(p))) return 0; 
        // One intersection point behind the trajectory starting point A:
        if (my_is_nan(intersection_s1(p))) return 1;
        // Both intersection points on the positive trajectory:
        return 2;
    }
#ifdef PRINTF_ENABLED
    printf("ERROR: THIS POINT SHOULD NOT BE REACHED. in number_of_intersections().\n");
#endif
}

floating_t squared_distance_from_center(floating_t X, floating_t Y, floating_t MX, floating_t MY)
{
    return (sqr(MX - X) + sqr(MY - Y));
}

bool intersecting_trajectory_starts_inside(IntersectionProblemParameters_t p)
{
    return (squared_distance_from_center(p.ax, p.ay, p.mx, p.my) < sqr(p.r));
}

inline bool intersecting_trajectory_starts_outside(IntersectionProblemParameters_t p)
{
    return ( ! intersecting_trajectory_starts_inside(p));
}

floating_t intersection_ratio_inside(IntersectionProblemParameters_t p)
{
    bool starts_inside = intersecting_trajectory_starts_inside(p);
    int num_of_intersections = number_of_intersections(p);
    
    // printf("HOLE ICE - INTERSECTION\n");
    // printf(" -> num of intersections = %i\n", num_of_intersections);
    // if (starts_inside) printf(" -> starts inside.\n");
    
    if (( ! starts_inside ) && ( num_of_intersections == 0 ))
        return 0.0;
    if (( ! starts_inside ) && ( num_of_intersections == 1 ))
        return 1.0 - intersection_s1(p);
    if (( ! starts_inside ) && ( num_of_intersections == 2 ))
        return intersection_s2(p) - intersection_s1(p);
    if (( starts_inside ) && ( num_of_intersections == 0 ))
        return 1.0;
    if (( starts_inside ) && ( num_of_intersections == 1 ))
        return intersection_s2(p);

#ifdef PRINTF_ENABLED
    printf("ERROR. This point should not be reached! in intersection_ratio_inside().\n");
#endif
    return my_nan();
}
