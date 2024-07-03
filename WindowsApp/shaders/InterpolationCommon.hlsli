
//===========================================================================
//
//                              INTERPOLATION
//
//===========================================================================

///Custom Fade In/Out through a point (x,y) in [0,1] range
float QuadraticThroughAGivenPoint(float x, float2 pointPos)
{
	float A = (1 - pointPos.y) / (1 - pointPos.x) - (pointPos.y / pointPos.x);
	float B = (A * (pointPos.x * pointPos.x) - pointPos.y) / pointPos.x;
	float y = A * (x * x) - B * (x);
	return y;
}

float3 LinearInterpolate(float t, float3 a, float3 b, float t0 = 0, float t1 = 1)
{
    //normalise to [0,1]
    t = (t - t0) / (t1 - t0);
    return a * (1 - t) + b * t;
}

float ConvertParameterToSmoothstep(float t)
{
    return 3 * pow(t, 2) - 2 * pow(t, 3);
}


float3 QuadraticInterpolate(
    float x, float3 y0, float3 y1, float3 y2,
    float x0 = 0, float x1 = 0.5f, float x2 = 1.f
 )
{
    return LinearInterpolate(x, LinearInterpolate(x, y0, y1, x0, x1), LinearInterpolate(x, y1, y2, x1, x2), x0, x2);
}

float3 QuadraticInterpolateBezier(float x, float3 y0, float3 y1, float3 y2, float x0 = 0, float x1 = 1.f)
{
    return LinearInterpolate(x, LinearInterpolate(x, y0, y1, x0, x1), LinearInterpolate(x, y1, y2, x0, x1), x0, x1);
}


float3 CubicInterpolate(float t, float3 a, float3 b, float3 c, float3 d, float t0 = 0, float t1 = 0.33f, float t2 = 0.66f, float t3 = 1.f)
{
    float3 quad1 = QuadraticInterpolate(t, a, b, c, t0, t1, t2);
    float3 quad2 = QuadraticInterpolate(t, b, c, d, t1, t2, t3);
    return LinearInterpolate(t, quad1, quad2, t0, t3);
}

float3 CubicInterpolateBezier(float t, float3 a, float3 b, float3 c, float3 d, float t0 = 0, float t1 = 1.f)
{
    float3 quad1 = QuadraticInterpolateBezier(t, a, b, c, t0, t1);
    float3 quad2 = QuadraticInterpolateBezier(t, b, c, d, t0, t1);
    return LinearInterpolate(t, quad1, quad2, t0, t1);
}

float3 InterpolateCubicSplineSegmentLagrange(float t, float3 sPrev, float3 s0, float3 s1, float3 sNext)
{
    float ti0 = -1.0;
    float ti1 =  0.0;
    float ti2 =  1.0;
    float ti3 =  2.0;
    //Lagrange Basis
    return
        sPrev * 
        (
            (t - ti1) / (ti0 - ti1) * 
            (t - ti2) / (ti0 - ti2) *
            (t - ti3) / (ti0 - ti3)
        ) +
        s0 * 
        (
            (t - ti0) / (ti1 - ti0) * 
            (t - ti2) / (ti1 - ti2) *
            (t - ti3) / (ti1 - ti3)
        ) +
        s1 * 
        (
            (t - ti0) / (ti2 - ti0) * 
            (t - ti1) / (ti2 - ti1) *
            (t - ti3) / (ti2 - ti3)
        ) +       
        sNext * 
        (
            (t - ti0) / (ti3 - ti0) * 
            (t - ti1) / (ti3 - ti1) *
            (t - ti2) / (ti3 - ti2)
        );
}

float3 InterpolateCubicSplineSegmentCatmullRom(float t, float3 A, float3 B, float3 C, float3 D)
{
    float t2 = t*t;
    float t3 = t*t*t;
    float3 a = -A/2.0 + (3.0*B)/2.0 - (3.0*C)/2.0 + D/2.0;
    float3 b = A - (5.0*B)/2.0 + 2.0*C - D / 2.0;
    float3 c = -A/2.0 + C/2.0;
    float3 d = B;
    
    return a*t3 + b*t2 + c*t + d;
}


float3 GetPointOnHermiteCurveParametric(float t, float3 y0, float3 d0, float3 y1, float3 d1)
{
    //Hermit Basis
    float H0 = 1 - 3 * pow(t, 2) + 2 * pow(t, 3);
    float H1 = t - 2 * pow(t, 2) + pow(t, 3);
    float H2 = 3 * pow(t, 2) - 2 * pow(t, 3);
    float H3 = -pow(t, 2) + pow(t, 3); 

    return (y0 * H0) + (d0 * H1) + (y1 * H2) + (d1 * H3);
}

float3 GetPointOnBezierCurveParametric(float t, float3 A, float3 B, float3 C, float3 D)
{
    //Bernstein Basis
    float B0 = pow(1 - t, 3);
    float B1 = 3 * t * pow(1 - t, 2);
    float B2 = 3 * t * t * (1 - t);
    float B3 = pow(t, 3);

    return A * B0 + B * B1 + C * B2 + D * B3;
}

float3 BilinearInterpolate(float2 xy, float3 s00, float3 s10, float3 s01, float3 s11, float t0 = 0, float t1 = 1)
{
    float3 horisontalPair1 = LinearInterpolate(xy.x, s00, s10);
    float3 horisontalPair2 = LinearInterpolate(xy.x, s01, s11);
    return LinearInterpolate(xy.y, horisontalPair1, horisontalPair2);
}



float MapToRange(float t, float t0, float t1, float newt0, float newt1)
{
	///Translate to origin, scale by ranges ratio, translate to new position
	return (t - t0) * ((newt1 - newt0) / (t1 - t0)) + newt0;
}

float2 MapToRange(float2 t, float t0, float t1, float newt0, float newt1)
{
	return float2(MapToRange(t.x, t0, t1, newt0, newt1), MapToRange(t.y, t0, t1, newt0, newt1));
}