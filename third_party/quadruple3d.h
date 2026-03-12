/* Automatically generated code, do not edit */
/* Generated from source file: quadruple3d.pck (DIM3 only) */

inline int quadruple_3d_filter( const double* p0, const double* p1, const double* p2, const double* p3) {
    double ac;
    ac = (((p0[0] * p2[0]) + (p0[1] * p2[1])) + (p0[2] * p2[2]));
    double bd;
    bd = (((p1[0] * p3[0]) + (p1[1] * p3[1])) + (p1[2] * p3[2]));
    double ad;
    ad = (((p0[0] * p3[0]) + (p0[1] * p3[1])) + (p0[2] * p3[2]));
    double bc;
    bc = (((p1[0] * p2[0]) + (p1[1] * p2[1])) + (p1[2] * p2[2]));
    double Delta;
    Delta = ((ac * bd) - (ad * bc));
    int int_tmp_result;
    double eps;
    double max1 = fabs(p0[0]);
    if( (max1 < fabs(p0[1])) )
    {
        max1 = fabs(p0[1]);
    } 
    if( (max1 < fabs(p0[2])) )
    {
        max1 = fabs(p0[2]);
    } 
    double max2 = fabs(p1[0]);
    if( (max2 < fabs(p1[1])) )
    {
        max2 = fabs(p1[1]);
    } 
    if( (max2 < fabs(p1[2])) )
    {
        max2 = fabs(p1[2]);
    } 
    double max3 = fabs(p2[0]);
    if( (max3 < fabs(p2[1])) )
    {
        max3 = fabs(p2[1]);
    } 
    if( (max3 < fabs(p2[2])) )
    {
        max3 = fabs(p2[2]);
    } 
    double max4 = fabs(p3[0]);
    if( (max4 < fabs(p3[1])) )
    {
        max4 = fabs(p3[1]);
    } 
    if( (max4 < fabs(p3[2])) )
    {
        max4 = fabs(p3[2]);
    } 
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    } 
    else 
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        } 
    } 
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    } 
    else 
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        } 
    } 
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    } 
    else 
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        } 
    } 
    if( (lower_bound_1 < 3.50952642631236329280e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    } 
    else 
    {
        if( (upper_bound_1 > 1.44740111546645180002e+76) )
        {
            return FPG_UNCERTAIN_VALUE;
        } 
        eps = (1.46673063709172807124e-14 * (((max1 * max3) * max2) * max4));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        } 
        else 
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            } 
            else 
            {
                return FPG_UNCERTAIN_VALUE;
            } 
        } 
    } 
    return int_tmp_result;
} 

