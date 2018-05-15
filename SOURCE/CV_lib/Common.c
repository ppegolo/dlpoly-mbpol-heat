double checkPBC (double x)
{
  if (x < 0.0)
    {
      do{
        x = x + 1;
      }while ( x < 0.0);
    }
  else if (x >= 1.0 )
    {
      do{
        x = x - 1;
      }while (x >= 1.0);
    }
  return x;
}

double checkPBC_Shifted (double x)
{
  if (x < -0.5)
    {
      do{
        x = x + 1;
      }while (x < -0.5);
    }
  else if (x > 0.5)
    {
      do{
        x = x - 1;
      }while (x > 0.5);
    }
  return x;
}


void GetInvH(double H[3][3], double invH[3][3])
{
  double detH;

  invH[0][0] = H[1][1]*H[2][2]-H[1][2]*H[2][1];
  invH[1][1] = H[2][2]*H[0][0]-H[2][0]*H[0][2];
  invH[2][2] = H[0][0]*H[1][1]-H[0][1]*H[1][0];

  invH[1][0] = H[1][2]*H[2][0]-H[1][0]*H[2][2];
  invH[2][1] = H[2][0]*H[0][1]-H[2][1]*H[0][0];
  invH[0][2] = H[0][1]*H[1][2]-H[0][2]*H[1][1];
  invH[2][0] = H[1][0]*H[2][1]-H[2][0]*H[1][1];
  invH[0][1] = H[2][1]*H[0][2]-H[0][1]*H[2][2];
  invH[1][2] = H[0][2]*H[1][0]-H[1][2]*H[0][0];

  detH = H[0][0]*invH[0][0] + H[0][1]*invH[1][0] + H[0][2]*invH[2][0];

  invH[0][0] /= detH;  invH[0][1] /= detH;  invH[0][2] /= detH;
  invH[1][0] /= detH;  invH[1][1] /= detH;  invH[1][2] /= detH;
  invH[2][0] /= detH;  invH[2][1] /= detH;  invH[2][2] /= detH;

  return ;
}


void Red2Cart(double x1[3], double H1[3][3], double x2[3])
{
  x2[0] = x1[0]*H1[0][0] + x1[1]*H1[1][0] + x1[2]*H1[2][0];
  x2[1] = x1[0]*H1[0][1] + x1[1]*H1[1][1] + x1[2]*H1[2][1];
  x2[2] = x1[0]*H1[0][2] + x1[1]*H1[1][2] + x1[2]*H1[2][2];

  return;
}

void Cart2Red(double x1[3], double invH1[3][3], double x2[3])
{
  x2[0] = x1[0]*invH1[0][0] + x1[1]*invH1[1][0] + x1[2]*invH1[2][0];
  x2[1] = x1[0]*invH1[0][1] + x1[1]*invH1[1][1] + x1[2]*invH1[2][1];
  x2[2] = x1[0]*invH1[0][2] + x1[1]*invH1[1][2] + x1[2]*invH1[2][2];

  return;
}

void V3V3subV3(double x1[3], double x2[3], double x3[3])
{
  x3[0] = x1[0] - x2[0];     
  x3[1] = x1[1] - x2[1];     
  x3[2] = x1[2] - x2[2];

  return;
}

double NormV3(double x1[3])
{
  double norm = x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2];
  return sqrt(norm);
}

void InitV3(double vec1[3])
{
  int i;
  for (i = 0; i < 3; i++)
    vec1[i] = 0.0;
  return;
}

void V3SmulV3(double x1[3], double x2, double x3[3])
{
  x3[0] = x1[0]*x2;   x3[1] = x1[1]*x2;   x3[2] = x1[2]*x2;
  return;
}












void InitV9(double vec1[9])
{
  int i;
  for (i = 0; i < 9; i++)
    vec1[i] = 0.0;
  return;
}

void InitV2(double vec1[2])
{
  vec1[0] = 0.00;
  vec1[1] = 0.00;
  return;
}

void V9symmetrizeV9(double vec1[9], double vec3[9])
{
  double vec2[9];
  int i;

  vec2[1] = 0.5*(vec1[1] + vec1[3]);
  vec2[3] = 0.5*(vec1[1] + vec1[3]);

  vec2[2] = 0.5*(vec1[6] + vec1[2]);
  vec2[6] = 0.5*(vec1[6] + vec1[2]);

  vec2[5] = 0.5*(vec1[5] + vec1[7]);
  vec2[7] = 0.5*(vec1[5] + vec1[7]);

  vec2[0] = vec1[0];
  vec2[4] = vec1[4];
  vec2[8] = vec1[8];

  for (i = 0; i < 9; i++)
    vec3[i] = vec2[i];

  return;
}

void V9symmetrizeV9addV9(double vec1[9], double vec3[9])
{
  double vec2[9];
  int i;

  vec2[1] = 0.5*(vec1[1] + vec1[3]);
  vec2[3] = 0.5*(vec1[1] + vec1[3]);

  vec2[2] = 0.5*(vec1[6] + vec1[2]);
  vec2[6] = 0.5*(vec1[6] + vec1[2]);

  vec2[5] = 0.5*(vec1[5] + vec1[7]);
  vec2[7] = 0.5*(vec1[5] + vec1[7]);

  vec2[0] = vec1[0];
  vec2[4] = vec1[4];
  vec2[8] = vec1[8];

  for (i = 0; i < 9; i++)
    vec3[i] += vec2[i];

  return;
}

void V3V3diadicV9(double vec1[3], double vec2[3], double vec3[9])
{
  vec3[0] = vec1[0]*vec2[0];
  vec3[1] = vec1[1]*vec2[0];
  vec3[2] = vec1[2]*vec2[0];

  vec3[3] = vec1[0]*vec2[1];
  vec3[4] = vec1[1]*vec2[1];
  vec3[5] = vec1[2]*vec2[1];

  vec3[6] = vec1[0]*vec2[2];
  vec3[7] = vec1[1]*vec2[2];
  vec3[8] = vec1[2]*vec2[2];
}

void V3V3diadicV9addV9(double vec1[3], double vec2[3], double vec3[9])
{
  vec3[0] += vec1[0]*vec2[0];
  vec3[1] += vec1[1]*vec2[0];
  vec3[2] += vec1[2]*vec2[0];

  vec3[3] += vec1[0]*vec2[1];
  vec3[4] += vec1[1]*vec2[1];
  vec3[5] += vec1[2]*vec2[1];

  vec3[6] += vec1[0]*vec2[2];
  vec3[7] += vec1[1]*vec2[2];
  vec3[8] += vec1[2]*vec2[2];
}

void V9V9addV9(double vec1[9], double vec2[9], double vec3[9])
{
  int i;
  for (i = 0; i < 9; i++)
    vec3[i] = vec1[i] + vec2[i];
}

void V9addV9(double vec1[9], double vec2[9])
{
  int i;
  for (i = 0; i < 9; i++)
    vec2[i] += vec1[i];
}

void Position_Difference_ReducedCoord(double x1[3], double x2[3], double x3[3])
{
  x3[0] = x1[0] - x2[0];     x3[1] = x1[1] - x2[1];     x3[2] = x1[2] - x2[2];
  return;
}

void V3M3mulV3(double x1[3], double H1[3][3], double x2[3])
{
  x2[0] = x1[0]*H1[0][0] + x1[1]*H1[1][0] + x1[2]*H1[2][0];
  x2[1] = x1[0]*H1[0][1] + x1[1]*H1[1][1] + x1[2]*H1[2][1];
  x2[2] = x1[0]*H1[0][2] + x1[1]*H1[1][2] + x1[2]*H1[2][2];

  return;
}

void Scalar_Multiply_V3(double x1[3], double x2, double x3[3])
{
  x3[0] = x1[0]*x2;   x3[1] = x1[1]*x2;   x3[2] = x1[2]*x2;
  return;
}

void SV3mulV3(double x1[3], double x2, double x3[3])
{
  x3[0] = x1[0]*x2;   x3[1] = x1[1]*x2;   x3[2] = x1[2]*x2;
  return;
}

void Scalar_Multiply_V3_Add(double x1[3], double x2, double x3[3])
{
  x3[0] += x1[0]*x2;   x3[1] += x1[1]*x2;   x3[2] += x1[2]*x2;
  return;
}

void V3SmulV3addV3(double x1[3], double x2, double x3[3])
{
  x3[0] += x1[0]*x2;   x3[1] += x1[1]*x2;   x3[2] += x1[2]*x2;
  return;
}

double NormSqV3(double x1[3])
{
  double norm = x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2];
  return norm;
}

double V3NormSqS(double x1[3])
{
  double norm = x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2];
  return norm;
}

void V3M9mulV3 (double *a1, double *a2, double *a3)
{ /* a1[3], a2[9], a3[3] */
  /* in a2[9] ordering is : 0 1 2; 3 4 5; 6 7 8;*/
  a3[0] = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2];
  a3[1] = a1[0]*a2[3] + a1[1]*a2[4] + a1[2]*a2[5];
  a3[2] = a1[0]*a2[6] + a1[1]*a2[7] + a1[2]*a2[8];
}

