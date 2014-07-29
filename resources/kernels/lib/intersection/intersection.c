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
    return (-intersection_beta(p) + 
        sign * my_sqrt(intersection_discriminant(p))) / 2 / intersection_alpha(p);
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
    else if (d > 0) return 2;
}


// int main()
// {
//     
//     // 1e3 Photonen, 1e3 Sprungpunkte, 1e2 Richtungen.
//     // Insgesamt dauert das 0,4 Sekunden. Prima.
//     for (int i = 0; i < 1e8; i++)
//     {
//     
//         IntersectionProblemParameters_t parameters = {
//             0.0, 0.5, // A
//             5.0, 0.5, // B
//             1.0, 1.0, // M
//             0.5
//         };
//         
//         // TODO: s-Parameterbereich auf [0;1] eingrenzen.
//         
//         // printf("Anzahl der Schnittpunkte: %i\n", number_of_intersections(parameters));
//         // printf("Schnittpunkt 1: (%e,%e)\n", X1(parameters), Y1(parameters));
//         // printf("Schnittpunkt 2: (%e,%e)\n", X2(parameters), Y2(parameters));
//         
//     }
//     return 0;
// }
