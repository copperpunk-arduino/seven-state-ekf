#pragma once
struct M66
{
	float m00 = 0.F;
	float m01 = 0.F;
	float m02 = 0.F;
	float m03 = 0.F;
	float m04 = 0.F;
	float m05 = 0.F;
	float m10 = 0.F;
	float m11 = 0.F;
	float m12 = 0.F;
	float m13 = 0.F;
	float m14 = 0.F;
	float m15 = 0.F;
	float m20 = 0.F;
	float m21 = 0.F;
	float m22 = 0.F;
	float m23 = 0.F;
	float m24 = 0.F;
	float m25 = 0.F;
	float m30 = 0.F;
	float m31 = 0.F;
	float m32 = 0.F;
	float m33 = 0.F;
	float m34 = 0.F;
	float m35 = 0.F;
	float m40 = 0.F;
	float m41 = 0.F;
	float m42 = 0.F;
	float m43 = 0.F;
	float m44 = 0.F;
	float m45 = 0.F;
	float m50 = 0.F;
	float m51 = 0.F;
	float m52 = 0.F;
	float m53 = 0.F;
	float m54 = 0.F;
	float m55 = 0.F;
};


class matrixSpec
{
public:
	static M66 invertMatrix(M66 m);
	static M66 invertMatrixFull(M66 m);
};

