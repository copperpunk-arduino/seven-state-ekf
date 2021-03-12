#include "MatrixSpec.h"

M66 matrixSpec::invertMatrix(M66 m) {
	float A4545 = m.m44 * m.m55 - m.m45 * m.m54;
	float A3545 = m.m43 * m.m55 - m.m45 * m.m53;
	float A3445 = m.m43 * m.m54 - m.m44 * m.m53;
	float A2545 = -m.m45 * m.m52;
	float A2445 = -m.m44 * m.m52;
	float A2345 = -m.m43 * m.m52;
	float A1545 = m.m41 * m.m55;
	float A1445 = m.m41 * m.m54;
	float A1345 = m.m41 * m.m53;
	float A1245 = m.m41 * m.m52;
	float A0545 = 0;
	float A0445 = 0;
	float A0345 = 0;
	float A0245 = 0;
	float A0145 = 0;
	float A4535 = m.m34 * m.m55 - m.m35 * m.m54;
	float A3535 = m.m33 * m.m55 - m.m35 * m.m53;
	float A3435 = m.m33 * m.m54 - m.m34 * m.m53;
	float A2535 = -m.m35 * m.m52;
	float A2435 = -m.m34 * m.m52;
	float A2335 = -m.m33 * m.m52;
	float A1535 = -m.m35 * m.m51;
	float A1435 = 0;
	float A1335 = 0;
	float A1235 = 0;
	float A4534 = m.m34 * m.m45 - m.m35 * m.m44;
	float A3534 = m.m33 * m.m45 - m.m35 * m.m43;
	float A3434 = m.m33 * m.m44 - m.m34 * m.m43;
	float A2534 = 0;
	float A2434 = 0;
	float A2334 = 0;
	float A1534 = -m.m35 * m.m41;
	float A1434 = -m.m34 * m.m41;
	float A1334 = -m.m33 * m.m41;
	float A1234 = 0;
	float A0535 = m.m30 * m.m55;
	float A0435 = m.m30 * m.m54;
	float A0335 = m.m30 * m.m53;
	float A0235 = m.m30 * m.m52;
	float A0534 = m.m30 * m.m45;
	float A0434 = m.m30 * m.m44;
	float A0334 = m.m30 * m.m43;
	float A0234 = 0;
	float A0135 = 0;
	float A0134 = m.m30 * m.m41;

	float B345345 = m.m33 * A4545 - m.m34 * A3545 + m.m35 * A3445;
	float B245345 = -m.m34 * A2545 + m.m35 * A2445;
	float B235345 = -m.m33 * A2545 + m.m35 * A2345;
	float B234345 = -m.m33 * A2445 + m.m34 * A2345;
	float B145345 = -m.m34 * A1545 + m.m35 * A1445;
	float B135345 = -m.m33 * A1545 + m.m35 * A1345;
	float B134345 = -m.m33 * A1445 + m.m34 * A1345;
	float B125345 = m.m35 * A1245;
	float B124345 = m.m34 * A1245;
	float B123345 = m.m33 * A1245;
	float B045345 = m.m30 * A4545 - m.m34 * A0545 + m.m35 * A0445;
	float B035345 = m.m30 * A3545 - m.m33 * A0545 + m.m35 * A0345;
	float B034345 = m.m30 * A3445 - m.m33 * A0445 + m.m34 * A0345;
	float B025345 = m.m30 * A2545 + m.m35 * A0245;
	float B024345 = m.m30 * A2445 + m.m34 * A0245;
	float B023345 = m.m30 * A2345 + m.m33 * A0245;
	float B015345 = m.m30 * A1545 + m.m35 * A0145;
	float B014345 = m.m30 * A1445 + m.m34 * A0145;
	float B013345 = m.m30 * A1345 + m.m33 * A0145;
	float B012345 = m.m30 * A1245;
	float B345245 = m.m25 * A3445;
	float B245245 = m.m22 * A4545 + m.m25 * A2445;
	float B235245 = m.m22 * A3545 + m.m25 * A2345;
	float B234245 = m.m22 * A3445;
	float B145245 = m.m25 * A1445;
	float B135245 = m.m25 * A1345;
	float B134245 = 0;
	float B125245 = -m.m22 * A1545 + m.m25 * A1245;
	float B124245 = -m.m22 * A1445;
	float B123245 = -m.m22 * A1345;
	float B345235 = m.m25 * A3435;
	float B245235 = m.m22 * A4535 + m.m25 * A2435;
	float B235235 = m.m22 * A3535 + m.m25 * A2335;
	float B234235 = m.m22 * A3435 + m.m24 * A2335;
	float B145235 = m.m25 * A1435;
	float B135235 = m.m25 * A1335;
	float B134235 = 0;
	float B125235 = -m.m22 * A1535 + m.m25 * A1235;
	float B124235 = -m.m22 * A1435;
	float B123235 = -m.m22 * A1335;
	float B345234 = m.m25 * A3434;
	float B245234 = m.m22 * A4534 + m.m25 * A2434;
	float B235234 = m.m22 * A3534 + m.m25 * A2334;
	float B234234 = m.m22 * A3434;
	float B145234 = m.m25 * A1434;
	float B135234 = m.m25 * A1334;
	float B134234 = 0;
	float B125234 = -m.m22 * A1534 + m.m25 * A1234;
	float B124234 = -m.m22 * A1434;
	float B123234 = -m.m22 * A1334;
	float B045245 = m.m25 * A0445;
	float B035245 = m.m25 * A0345;
	float B034245 = 0;
	float B025245 = -m.m22 * A0545 + m.m25 * A0245;
	float B024245 = -m.m22 * A0445;
	float B023245 = -m.m22 * A0345;
	float B045235 = m.m25 * A0435;
	float B035235 = m.m25 * A0335;
	float B034235 = 0;
	float B025235 = -m.m22 * A0535 + m.m25 * A0235;
	float B024235 = -m.m22 * A0435;
	float B023235 = -m.m22 * A0335;
	float B045234 = m.m25 * A0434;
	float B035234 = m.m25 * A0334;
	float B034234 = 0;
	float B025234 = -m.m22 * A0534 + m.m25 * A0234;
	float B024234 = -m.m22 * A0434;
	float B023234 = -m.m22 * A0334;
	float B015245 = m.m25 * A0145;
	float B014245 = 0;
	float B013245 = 0;
	float B015235 = m.m25 * A0135;
	float B014235 = 0;
	float B013235 = 0;
	float B015234 = m.m25 * A0134;
	float B014234 = 0;
	float B013234 = 0;
	float B012245 = m.m22 * A0145;
	float B012235 = m.m22 * A0135;
	float B012234 = m.m22 * A0134;

	float C23452345 = m.m22 * B345345 - m.m25 * B234345;
	float C13452345 = -m.m25 * B134345;
	float C12452345 = -m.m22 * B145345 - m.m25 * B124345;
	float C12352345 = -m.m22 * B135345 - m.m25 * B123345;
	float C12342345 = -m.m22 * B134345;
	float C03452345 = -m.m25 * B034345;
	float C02452345 = -m.m22 * B045345 - m.m25 * B024345;
	float C02352345 = -m.m22 * B035345 - m.m25 * B023345;
	float C02342345 = -m.m22 * B034345;
	float C01452345 = -m.m25 * B014345;
	float C01352345 = -m.m25 * B013345;
	float C01252345 = m.m22 * B015345 - m.m25 * B012345;
	float C01242345 = m.m22 * B014345;
	float C01232345 = m.m22 * B013345;
	float C23451345 = m.m14 * B235345;
	float C13451345 = m.m11 * B345345 + m.m14 * B135345;
	float C12451345 = m.m11 * B245345 + m.m14 * B125345;
	float C12351345 = m.m11 * B235345;
	float C12341345 = m.m11 * B234345 - m.m14 * B123345;
	float C23451245 = m.m14 * B235245;
	float C13451245 = m.m11 * B345245 + m.m14 * B135245;
	float C12451245 = m.m11 * B245245 + m.m14 * B125245;
	float C12351245 = m.m11 * B235245;
	float C12341245 = m.m11 * B234245 - m.m14 * B123245;
	float C23451235 = m.m14 * B235235;
	float C13451235 = m.m11 * B345235 + m.m14 * B135235;
	float C12451235 = m.m11 * B245235 + m.m14 * B125235;
	float C12351235 = m.m11 * B235235;
	float C12341235 = m.m11 * B234235 - m.m14 * B123235;
	float C23451234 = m.m14 * B235234;
	float C13451234 = m.m11 * B345234 + m.m14 * B135234;
	float C12451234 = m.m11 * B245234 + m.m14 * B125234;
	float C12351234 = m.m11 * B235234;
	float C12341234 = m.m11 * B234234 - m.m14 * B123234;
	float C03451345 = m.m14 * B035345;
	float C02451345 = m.m14 * B025345;
	float C03451245 = m.m14 * B035245;
	float C02451245 = m.m14 * B025245;
	float C02341245 = -m.m14 * B023245;
	float C03451235 = m.m14 * B035235;
	float C02451235 = m.m14 * B025235;
	float C02341235 = -m.m14 * B023235;
	float C03451234 = m.m14 * B035234;
	float C02451234 = m.m14 * B025234;
	float C02341234 = -m.m14 * B023234;
	float C01451345 = -m.m11 * B045345 + m.m14 * B015345;
	float C01351345 = -m.m11 * B035345;
	float C01341345 = -m.m11 * B034345 - m.m14 * B013345;
	float C01451245 = -m.m11 * B045245 + m.m14 * B015245;
	float C01351245 = -m.m11 * B035245;
	float C01341245 = -m.m11 * B034245 - m.m14 * B013245;
	float C01451235 = -m.m11 * B045235 + m.m14 * B015235;
	float C01351235 = -m.m11 * B035235;
	float C01341235 = -m.m11 * B034235 - m.m14 * B013235;
	float C01451234 = -m.m11 * B045234 + m.m14 * B015234;
	float C01351234 = -m.m11 * B035234;
	float C01341234 = -m.m11 * B034234 - m.m14 * B013234;
	float C01251345 = -m.m11 * B025345;
	float C01241345 = -m.m11 * B024345 - m.m14 * B012345;
	float C01251245 = -m.m11 * B025245;
	float C01241245 = -m.m11 * B024245 - m.m14 * B012245;
	float C01251235 = -m.m11 * B025235;
	float C01241235 = -m.m11 * B024235 - m.m14 * B012235;
	float C01251234 = -m.m11 * B025234;
	float C01241234 = -m.m11 * B024234 - m.m14 * B012234;
	float C01231345 = -m.m11 * B023345;
	float C01231245 = -m.m11 * B023245;
	float C01231235 = -m.m11 * B023235;
	float C01231234 = -m.m11 * B023234;

	float det = m.m00 * (m.m11 * C23452345 - m.m14 * C12352345)
		- m.m03 * (-m.m11 * C02452345 - m.m14 * C01252345);

	if (det != 0.F) {
		det = 1.0F / det;
	}

	M66 newM;
	newM.m00 = det * (m.m11 * C23452345 - m.m14 * C12352345);
	newM.m01 = det * -(m.m03 * C12452345);
	newM.m02 = det * (m.m03 * C12451345);
	newM.m03 = det * -(m.m03 * C12451245);
	newM.m04 = det * (m.m03 * C12451235);
	newM.m05 = det * -(m.m03 * C12451234);
	newM.m10 = det * -(-m.m14 * C02352345);
	newM.m11 = det * (m.m00 * C23452345 + m.m03 * C02452345);
	newM.m12 = det * -(m.m00 * C23451345 + m.m03 * C02451345);
	newM.m13 = det * (m.m00 * C23451245 + m.m03 * C02451245);
	newM.m14 = det * -(m.m00 * C23451235 + m.m03 * C02451235);
	newM.m15 = det * (m.m00 * C23451234 + m.m03 * C02451234);
	newM.m20 = det * (-m.m11 * C03452345 - m.m14 * C01352345);
	newM.m21 = det * -(m.m00 * C13452345 + m.m03 * C01452345);
	newM.m22 = det * (m.m00 * C13451345 + m.m03 * C01451345);
	newM.m23 = det * -(m.m00 * C13451245 + m.m03 * C01451245);
	newM.m24 = det * (m.m00 * C13451235 + m.m03 * C01451235);
	newM.m25 = det * -(m.m00 * C13451234 + m.m03 * C01451234);
	newM.m30 = det * -(-m.m11 * C02452345 - m.m14 * C01252345);
	newM.m31 = det * (m.m00 * C12452345);
	newM.m32 = det * -(m.m00 * C12451345);
	newM.m33 = det * (m.m00 * C12451245);
	newM.m34 = det * -(m.m00 * C12451235);
	newM.m35 = det * (m.m00 * C12451234);
	newM.m40 = det * (-m.m11 * C02352345);
	newM.m41 = det * -(m.m00 * C12352345 - m.m03 * C01252345);
	newM.m42 = det * (m.m00 * C12351345 - m.m03 * C01251345);
	newM.m43 = det * -(m.m00 * C12351245 - m.m03 * C01251245);
	newM.m44 = det * (m.m00 * C12351235 - m.m03 * C01251235);
	newM.m45 = det * -(m.m00 * C12351234 - m.m03 * C01251234);
	newM.m50 = det * -(-m.m11 * C02342345 + m.m14 * C01232345);
	newM.m51 = det * (m.m00 * C12342345 - m.m03 * C01242345);
	newM.m52 = det * -(m.m00 * C12341345 - m.m03 * C01241345);
	newM.m53 = det * (m.m00 * C12341245 - m.m03 * C01241245);
	newM.m54 = det * -(m.m00 * C12341235 - m.m03 * C01241235);
	newM.m55 = det * (m.m00 * C12341234 - m.m03 * C01241234);
	return newM;
}

M66 matrixSpec::invertMatrixFull(M66 m) {

	float A4545 = m.m44 * m.m55 - m.m45 * m.m54;
	float A3545 = m.m43 * m.m55 - m.m45 * m.m53;
	float A3445 = m.m43 * m.m54 - m.m44 * m.m53;
	float A2545 = m.m42 * m.m55 - m.m45 * m.m52;
	float A2445 = m.m42 * m.m54 - m.m44 * m.m52;
	float A2345 = m.m42 * m.m53 - m.m43 * m.m52;
	float A1545 = m.m41 * m.m55 - m.m45 * m.m51;
	float A1445 = m.m41 * m.m54 - m.m44 * m.m51;
	float A1345 = m.m41 * m.m53 - m.m43 * m.m51;
	float A1245 = m.m41 * m.m52 - m.m42 * m.m51;
	float A0545 = m.m40 * m.m55 - m.m45 * m.m50;
	float A0445 = m.m40 * m.m54 - m.m44 * m.m50;
	float A0345 = m.m40 * m.m53 - m.m43 * m.m50;
	float A0245 = m.m40 * m.m52 - m.m42 * m.m50;
	float A0145 = m.m40 * m.m51 - m.m41 * m.m50;
	float A4535 = m.m34 * m.m55 - m.m35 * m.m54;
	float A3535 = m.m33 * m.m55 - m.m35 * m.m53;
	float A3435 = m.m33 * m.m54 - m.m34 * m.m53;
	float A2535 = m.m32 * m.m55 - m.m35 * m.m52;
	float A2435 = m.m32 * m.m54 - m.m34 * m.m52;
	float A2335 = m.m32 * m.m53 - m.m33 * m.m52;
	float A1535 = m.m31 * m.m55 - m.m35 * m.m51;
	float A1435 = m.m31 * m.m54 - m.m34 * m.m51;
	float A1335 = m.m31 * m.m53 - m.m33 * m.m51;
	float A1235 = m.m31 * m.m52 - m.m32 * m.m51;
	float A4534 = m.m34 * m.m45 - m.m35 * m.m44;
	float A3534 = m.m33 * m.m45 - m.m35 * m.m43;
	float A3434 = m.m33 * m.m44 - m.m34 * m.m43;
	float A2534 = m.m32 * m.m45 - m.m35 * m.m42;
	float A2434 = m.m32 * m.m44 - m.m34 * m.m42;
	float A2334 = m.m32 * m.m43 - m.m33 * m.m42;
	float A1534 = m.m31 * m.m45 - m.m35 * m.m41;
	float A1434 = m.m31 * m.m44 - m.m34 * m.m41;
	float A1334 = m.m31 * m.m43 - m.m33 * m.m41;
	float A1234 = m.m31 * m.m42 - m.m32 * m.m41;
	float A0535 = m.m30 * m.m55 - m.m35 * m.m50;
	float A0435 = m.m30 * m.m54 - m.m34 * m.m50;
	float A0335 = m.m30 * m.m53 - m.m33 * m.m50;
	float A0235 = m.m30 * m.m52 - m.m32 * m.m50;
	float A0534 = m.m30 * m.m45 - m.m35 * m.m40;
	float A0434 = m.m30 * m.m44 - m.m34 * m.m40;
	float A0334 = m.m30 * m.m43 - m.m33 * m.m40;
	float A0234 = m.m30 * m.m42 - m.m32 * m.m40;
	float A0135 = m.m30 * m.m51 - m.m31 * m.m50;
	float A0134 = m.m30 * m.m41 - m.m31 * m.m40;

	float B345345 = m.m33 * A4545 - m.m34 * A3545 + m.m35 * A3445;
	float B245345 = m.m32 * A4545 - m.m34 * A2545 + m.m35 * A2445;
	float B235345 = m.m32 * A3545 - m.m33 * A2545 + m.m35 * A2345;
	float B234345 = m.m32 * A3445 - m.m33 * A2445 + m.m34 * A2345;
	float B145345 = m.m31 * A4545 - m.m34 * A1545 + m.m35 * A1445;
	float B135345 = m.m31 * A3545 - m.m33 * A1545 + m.m35 * A1345;
	float B134345 = m.m31 * A3445 - m.m33 * A1445 + m.m34 * A1345;
	float B125345 = m.m31 * A2545 - m.m32 * A1545 + m.m35 * A1245;
	float B124345 = m.m31 * A2445 - m.m32 * A1445 + m.m34 * A1245;
	float B123345 = m.m31 * A2345 - m.m32 * A1345 + m.m33 * A1245;
	float B045345 = m.m30 * A4545 - m.m34 * A0545 + m.m35 * A0445;
	float B035345 = m.m30 * A3545 - m.m33 * A0545 + m.m35 * A0345;
	float B034345 = m.m30 * A3445 - m.m33 * A0445 + m.m34 * A0345;
	float B025345 = m.m30 * A2545 - m.m32 * A0545 + m.m35 * A0245;
	float B024345 = m.m30 * A2445 - m.m32 * A0445 + m.m34 * A0245;
	float B023345 = m.m30 * A2345 - m.m32 * A0345 + m.m33 * A0245;
	float B015345 = m.m30 * A1545 - m.m31 * A0545 + m.m35 * A0145;
	float B014345 = m.m30 * A1445 - m.m31 * A0445 + m.m34 * A0145;
	float B013345 = m.m30 * A1345 - m.m31 * A0345 + m.m33 * A0145;
	float B012345 = m.m30 * A1245 - m.m31 * A0245 + m.m32 * A0145;
	float B345245 = m.m23 * A4545 - m.m24 * A3545 + m.m25 * A3445;
	float B245245 = m.m22 * A4545 - m.m24 * A2545 + m.m25 * A2445;
	float B235245 = m.m22 * A3545 - m.m23 * A2545 + m.m25 * A2345;
	float B234245 = m.m22 * A3445 - m.m23 * A2445 + m.m24 * A2345;
	float B145245 = m.m21 * A4545 - m.m24 * A1545 + m.m25 * A1445;
	float B135245 = m.m21 * A3545 - m.m23 * A1545 + m.m25 * A1345;
	float B134245 = m.m21 * A3445 - m.m23 * A1445 + m.m24 * A1345;
	float B125245 = m.m21 * A2545 - m.m22 * A1545 + m.m25 * A1245;
	float B124245 = m.m21 * A2445 - m.m22 * A1445 + m.m24 * A1245;
	float B123245 = m.m21 * A2345 - m.m22 * A1345 + m.m23 * A1245;
	float B345235 = m.m23 * A4535 - m.m24 * A3535 + m.m25 * A3435;
	float B245235 = m.m22 * A4535 - m.m24 * A2535 + m.m25 * A2435;
	float B235235 = m.m22 * A3535 - m.m23 * A2535 + m.m25 * A2335;
	float B234235 = m.m22 * A3435 - m.m23 * A2435 + m.m24 * A2335;
	float B145235 = m.m21 * A4535 - m.m24 * A1535 + m.m25 * A1435;
	float B135235 = m.m21 * A3535 - m.m23 * A1535 + m.m25 * A1335;
	float B134235 = m.m21 * A3435 - m.m23 * A1435 + m.m24 * A1335;
	float B125235 = m.m21 * A2535 - m.m22 * A1535 + m.m25 * A1235;
	float B124235 = m.m21 * A2435 - m.m22 * A1435 + m.m24 * A1235;
	float B123235 = m.m21 * A2335 - m.m22 * A1335 + m.m23 * A1235;
	float B345234 = m.m23 * A4534 - m.m24 * A3534 + m.m25 * A3434;
	float B245234 = m.m22 * A4534 - m.m24 * A2534 + m.m25 * A2434;
	float B235234 = m.m22 * A3534 - m.m23 * A2534 + m.m25 * A2334;
	float B234234 = m.m22 * A3434 - m.m23 * A2434 + m.m24 * A2334;
	float B145234 = m.m21 * A4534 - m.m24 * A1534 + m.m25 * A1434;
	float B135234 = m.m21 * A3534 - m.m23 * A1534 + m.m25 * A1334;
	float B134234 = m.m21 * A3434 - m.m23 * A1434 + m.m24 * A1334;
	float B125234 = m.m21 * A2534 - m.m22 * A1534 + m.m25 * A1234;
	float B124234 = m.m21 * A2434 - m.m22 * A1434 + m.m24 * A1234;
	float B123234 = m.m21 * A2334 - m.m22 * A1334 + m.m23 * A1234;
	float B045245 = m.m20 * A4545 - m.m24 * A0545 + m.m25 * A0445;
	float B035245 = m.m20 * A3545 - m.m23 * A0545 + m.m25 * A0345;
	float B034245 = m.m20 * A3445 - m.m23 * A0445 + m.m24 * A0345;
	float B025245 = m.m20 * A2545 - m.m22 * A0545 + m.m25 * A0245;
	float B024245 = m.m20 * A2445 - m.m22 * A0445 + m.m24 * A0245;
	float B023245 = m.m20 * A2345 - m.m22 * A0345 + m.m23 * A0245;
	float B045235 = m.m20 * A4535 - m.m24 * A0535 + m.m25 * A0435;
	float B035235 = m.m20 * A3535 - m.m23 * A0535 + m.m25 * A0335;
	float B034235 = m.m20 * A3435 - m.m23 * A0435 + m.m24 * A0335;
	float B025235 = m.m20 * A2535 - m.m22 * A0535 + m.m25 * A0235;
	float B024235 = m.m20 * A2435 - m.m22 * A0435 + m.m24 * A0235;
	float B023235 = m.m20 * A2335 - m.m22 * A0335 + m.m23 * A0235;
	float B045234 = m.m20 * A4534 - m.m24 * A0534 + m.m25 * A0434;
	float B035234 = m.m20 * A3534 - m.m23 * A0534 + m.m25 * A0334;
	float B034234 = m.m20 * A3434 - m.m23 * A0434 + m.m24 * A0334;
	float B025234 = m.m20 * A2534 - m.m22 * A0534 + m.m25 * A0234;
	float B024234 = m.m20 * A2434 - m.m22 * A0434 + m.m24 * A0234;
	float B023234 = m.m20 * A2334 - m.m22 * A0334 + m.m23 * A0234;
	float B015245 = m.m20 * A1545 - m.m21 * A0545 + m.m25 * A0145;
	float B014245 = m.m20 * A1445 - m.m21 * A0445 + m.m24 * A0145;
	float B013245 = m.m20 * A1345 - m.m21 * A0345 + m.m23 * A0145;
	float B015235 = m.m20 * A1535 - m.m21 * A0535 + m.m25 * A0135;
	float B014235 = m.m20 * A1435 - m.m21 * A0435 + m.m24 * A0135;
	float B013235 = m.m20 * A1335 - m.m21 * A0335 + m.m23 * A0135;
	float B015234 = m.m20 * A1534 - m.m21 * A0534 + m.m25 * A0134;
	float B014234 = m.m20 * A1434 - m.m21 * A0434 + m.m24 * A0134;
	float B013234 = m.m20 * A1334 - m.m21 * A0334 + m.m23 * A0134;
	float B012245 = m.m20 * A1245 - m.m21 * A0245 + m.m22 * A0145;
	float B012235 = m.m20 * A1235 - m.m21 * A0235 + m.m22 * A0135;
	float B012234 = m.m20 * A1234 - m.m21 * A0234 + m.m22 * A0134;

	float C23452345 = m.m22 * B345345 - m.m23 * B245345 + m.m24 * B235345 - m.m25 * B234345;
	float C13452345 = m.m21 * B345345 - m.m23 * B145345 + m.m24 * B135345 - m.m25 * B134345;
	float C12452345 = m.m21 * B245345 - m.m22 * B145345 + m.m24 * B125345 - m.m25 * B124345;
	float C12352345 = m.m21 * B235345 - m.m22 * B135345 + m.m23 * B125345 - m.m25 * B123345;
	float C12342345 = m.m21 * B234345 - m.m22 * B134345 + m.m23 * B124345 - m.m24 * B123345;
	float C03452345 = m.m20 * B345345 - m.m23 * B045345 + m.m24 * B035345 - m.m25 * B034345;
	float C02452345 = m.m20 * B245345 - m.m22 * B045345 + m.m24 * B025345 - m.m25 * B024345;
	float C02352345 = m.m20 * B235345 - m.m22 * B035345 + m.m23 * B025345 - m.m25 * B023345;
	float C02342345 = m.m20 * B234345 - m.m22 * B034345 + m.m23 * B024345 - m.m24 * B023345;
	float C01452345 = m.m20 * B145345 - m.m21 * B045345 + m.m24 * B015345 - m.m25 * B014345;
	float C01352345 = m.m20 * B135345 - m.m21 * B035345 + m.m23 * B015345 - m.m25 * B013345;
	float C01342345 = m.m20 * B134345 - m.m21 * B034345 + m.m23 * B014345 - m.m24 * B013345;
	float C01252345 = m.m20 * B125345 - m.m21 * B025345 + m.m22 * B015345 - m.m25 * B012345;
	float C01242345 = m.m20 * B124345 - m.m21 * B024345 + m.m22 * B014345 - m.m24 * B012345;
	float C01232345 = m.m20 * B123345 - m.m21 * B023345 + m.m22 * B013345 - m.m23 * B012345;
	float C23451345 = m.m12 * B345345 - m.m13 * B245345 + m.m14 * B235345 - m.m15 * B234345;
	float C13451345 = m.m11 * B345345 - m.m13 * B145345 + m.m14 * B135345 - m.m15 * B134345;
	float C12451345 = m.m11 * B245345 - m.m12 * B145345 + m.m14 * B125345 - m.m15 * B124345;
	float C12351345 = m.m11 * B235345 - m.m12 * B135345 + m.m13 * B125345 - m.m15 * B123345;
	float C12341345 = m.m11 * B234345 - m.m12 * B134345 + m.m13 * B124345 - m.m14 * B123345;
	float C23451245 = m.m12 * B345245 - m.m13 * B245245 + m.m14 * B235245 - m.m15 * B234245;
	float C13451245 = m.m11 * B345245 - m.m13 * B145245 + m.m14 * B135245 - m.m15 * B134245;
	float C12451245 = m.m11 * B245245 - m.m12 * B145245 + m.m14 * B125245 - m.m15 * B124245;
	float C12351245 = m.m11 * B235245 - m.m12 * B135245 + m.m13 * B125245 - m.m15 * B123245;
	float C12341245 = m.m11 * B234245 - m.m12 * B134245 + m.m13 * B124245 - m.m14 * B123245;
	float C23451235 = m.m12 * B345235 - m.m13 * B245235 + m.m14 * B235235 - m.m15 * B234235;
	float C13451235 = m.m11 * B345235 - m.m13 * B145235 + m.m14 * B135235 - m.m15 * B134235;
	float C12451235 = m.m11 * B245235 - m.m12 * B145235 + m.m14 * B125235 - m.m15 * B124235;
	float C12351235 = m.m11 * B235235 - m.m12 * B135235 + m.m13 * B125235 - m.m15 * B123235;
	float C12341235 = m.m11 * B234235 - m.m12 * B134235 + m.m13 * B124235 - m.m14 * B123235;
	float C23451234 = m.m12 * B345234 - m.m13 * B245234 + m.m14 * B235234 - m.m15 * B234234;
	float C13451234 = m.m11 * B345234 - m.m13 * B145234 + m.m14 * B135234 - m.m15 * B134234;
	float C12451234 = m.m11 * B245234 - m.m12 * B145234 + m.m14 * B125234 - m.m15 * B124234;
	float C12351234 = m.m11 * B235234 - m.m12 * B135234 + m.m13 * B125234 - m.m15 * B123234;
	float C12341234 = m.m11 * B234234 - m.m12 * B134234 + m.m13 * B124234 - m.m14 * B123234;
	float C03451345 = m.m10 * B345345 - m.m13 * B045345 + m.m14 * B035345 - m.m15 * B034345;
	float C02451345 = m.m10 * B245345 - m.m12 * B045345 + m.m14 * B025345 - m.m15 * B024345;
	float C02351345 = m.m10 * B235345 - m.m12 * B035345 + m.m13 * B025345 - m.m15 * B023345;
	float C02341345 = m.m10 * B234345 - m.m12 * B034345 + m.m13 * B024345 - m.m14 * B023345;
	float C03451245 = m.m10 * B345245 - m.m13 * B045245 + m.m14 * B035245 - m.m15 * B034245;
	float C02451245 = m.m10 * B245245 - m.m12 * B045245 + m.m14 * B025245 - m.m15 * B024245;
	float C02351245 = m.m10 * B235245 - m.m12 * B035245 + m.m13 * B025245 - m.m15 * B023245;
	float C02341245 = m.m10 * B234245 - m.m12 * B034245 + m.m13 * B024245 - m.m14 * B023245;
	float C03451235 = m.m10 * B345235 - m.m13 * B045235 + m.m14 * B035235 - m.m15 * B034235;
	float C02451235 = m.m10 * B245235 - m.m12 * B045235 + m.m14 * B025235 - m.m15 * B024235;
	float C02351235 = m.m10 * B235235 - m.m12 * B035235 + m.m13 * B025235 - m.m15 * B023235;
	float C02341235 = m.m10 * B234235 - m.m12 * B034235 + m.m13 * B024235 - m.m14 * B023235;
	float C03451234 = m.m10 * B345234 - m.m13 * B045234 + m.m14 * B035234 - m.m15 * B034234;
	float C02451234 = m.m10 * B245234 - m.m12 * B045234 + m.m14 * B025234 - m.m15 * B024234;
	float C02351234 = m.m10 * B235234 - m.m12 * B035234 + m.m13 * B025234 - m.m15 * B023234;
	float C02341234 = m.m10 * B234234 - m.m12 * B034234 + m.m13 * B024234 - m.m14 * B023234;
	float C01451345 = m.m10 * B145345 - m.m11 * B045345 + m.m14 * B015345 - m.m15 * B014345;
	float C01351345 = m.m10 * B135345 - m.m11 * B035345 + m.m13 * B015345 - m.m15 * B013345;
	float C01341345 = m.m10 * B134345 - m.m11 * B034345 + m.m13 * B014345 - m.m14 * B013345;
	float C01451245 = m.m10 * B145245 - m.m11 * B045245 + m.m14 * B015245 - m.m15 * B014245;
	float C01351245 = m.m10 * B135245 - m.m11 * B035245 + m.m13 * B015245 - m.m15 * B013245;
	float C01341245 = m.m10 * B134245 - m.m11 * B034245 + m.m13 * B014245 - m.m14 * B013245;
	float C01451235 = m.m10 * B145235 - m.m11 * B045235 + m.m14 * B015235 - m.m15 * B014235;
	float C01351235 = m.m10 * B135235 - m.m11 * B035235 + m.m13 * B015235 - m.m15 * B013235;
	float C01341235 = m.m10 * B134235 - m.m11 * B034235 + m.m13 * B014235 - m.m14 * B013235;
	float C01451234 = m.m10 * B145234 - m.m11 * B045234 + m.m14 * B015234 - m.m15 * B014234;
	float C01351234 = m.m10 * B135234 - m.m11 * B035234 + m.m13 * B015234 - m.m15 * B013234;
	float C01341234 = m.m10 * B134234 - m.m11 * B034234 + m.m13 * B014234 - m.m14 * B013234;
	float C01251345 = m.m10 * B125345 - m.m11 * B025345 + m.m12 * B015345 - m.m15 * B012345;
	float C01241345 = m.m10 * B124345 - m.m11 * B024345 + m.m12 * B014345 - m.m14 * B012345;
	float C01251245 = m.m10 * B125245 - m.m11 * B025245 + m.m12 * B015245 - m.m15 * B012245;
	float C01241245 = m.m10 * B124245 - m.m11 * B024245 + m.m12 * B014245 - m.m14 * B012245;
	float C01251235 = m.m10 * B125235 - m.m11 * B025235 + m.m12 * B015235 - m.m15 * B012235;
	float C01241235 = m.m10 * B124235 - m.m11 * B024235 + m.m12 * B014235 - m.m14 * B012235;
	float C01251234 = m.m10 * B125234 - m.m11 * B025234 + m.m12 * B015234 - m.m15 * B012234;
	float C01241234 = m.m10 * B124234 - m.m11 * B024234 + m.m12 * B014234 - m.m14 * B012234;
	float C01231345 = m.m10 * B123345 - m.m11 * B023345 + m.m12 * B013345 - m.m13 * B012345;
	float C01231245 = m.m10 * B123245 - m.m11 * B023245 + m.m12 * B013245 - m.m13 * B012245;
	float C01231235 = m.m10 * B123235 - m.m11 * B023235 + m.m12 * B013235 - m.m13 * B012235;
	float C01231234 = m.m10 * B123234 - m.m11 * B023234 + m.m12 * B013234 - m.m13 * B012234;

	float det = m.m00 * (m.m11 * C23452345 - m.m12 * C13452345 + m.m13 * C12452345 - m.m14 * C12352345 + m.m15 * C12342345)
		- m.m01 * (m.m10 * C23452345 - m.m12 * C03452345 + m.m13 * C02452345 - m.m14 * C02352345 + m.m15 * C02342345)
		+ m.m02 * (m.m10 * C13452345 - m.m11 * C03452345 + m.m13 * C01452345 - m.m14 * C01352345 + m.m15 * C01342345)
		- m.m03 * (m.m10 * C12452345 - m.m11 * C02452345 + m.m12 * C01452345 - m.m14 * C01252345 + m.m15 * C01242345)
		+ m.m04 * (m.m10 * C12352345 - m.m11 * C02352345 + m.m12 * C01352345 - m.m13 * C01252345 + m.m15 * C01232345)
		- m.m05 * (m.m10 * C12342345 - m.m11 * C02342345 + m.m12 * C01342345 - m.m13 * C01242345 + m.m14 * C01232345);
	if (det != 0) {
		det = 1 / det;
	}

	M66 newM;

	newM.m00 = det * (m.m11 * C23452345 - m.m12 * C13452345 + m.m13 * C12452345 - m.m14 * C12352345 + m.m15 * C12342345);
	newM.m01 = det * -(m.m01 * C23452345 - m.m02 * C13452345 + m.m03 * C12452345 - m.m04 * C12352345 + m.m05 * C12342345);
	newM.m02 = det * (m.m01 * C23451345 - m.m02 * C13451345 + m.m03 * C12451345 - m.m04 * C12351345 + m.m05 * C12341345);
	newM.m03 = det * -(m.m01 * C23451245 - m.m02 * C13451245 + m.m03 * C12451245 - m.m04 * C12351245 + m.m05 * C12341245);
	newM.m04 = det * (m.m01 * C23451235 - m.m02 * C13451235 + m.m03 * C12451235 - m.m04 * C12351235 + m.m05 * C12341235);
	newM.m05 = det * -(m.m01 * C23451234 - m.m02 * C13451234 + m.m03 * C12451234 - m.m04 * C12351234 + m.m05 * C12341234);
	newM.m10 = det * -(m.m10 * C23452345 - m.m12 * C03452345 + m.m13 * C02452345 - m.m14 * C02352345 + m.m15 * C02342345);
	newM.m11 = det * (m.m00 * C23452345 - m.m02 * C03452345 + m.m03 * C02452345 - m.m04 * C02352345 + m.m05 * C02342345);
	newM.m12 = det * -(m.m00 * C23451345 - m.m02 * C03451345 + m.m03 * C02451345 - m.m04 * C02351345 + m.m05 * C02341345);
	newM.m13 = det * (m.m00 * C23451245 - m.m02 * C03451245 + m.m03 * C02451245 - m.m04 * C02351245 + m.m05 * C02341245);
	newM.m14 = det * -(m.m00 * C23451235 - m.m02 * C03451235 + m.m03 * C02451235 - m.m04 * C02351235 + m.m05 * C02341235);
	newM.m15 = det * (m.m00 * C23451234 - m.m02 * C03451234 + m.m03 * C02451234 - m.m04 * C02351234 + m.m05 * C02341234);
	newM.m20 = det * (m.m10 * C13452345 - m.m11 * C03452345 + m.m13 * C01452345 - m.m14 * C01352345 + m.m15 * C01342345);
	newM.m21 = det * -(m.m00 * C13452345 - m.m01 * C03452345 + m.m03 * C01452345 - m.m04 * C01352345 + m.m05 * C01342345);
	newM.m22 = det * (m.m00 * C13451345 - m.m01 * C03451345 + m.m03 * C01451345 - m.m04 * C01351345 + m.m05 * C01341345);
	newM.m23 = det * -(m.m00 * C13451245 - m.m01 * C03451245 + m.m03 * C01451245 - m.m04 * C01351245 + m.m05 * C01341245);
	newM.m24 = det * (m.m00 * C13451235 - m.m01 * C03451235 + m.m03 * C01451235 - m.m04 * C01351235 + m.m05 * C01341235);
	newM.m25 = det * -(m.m00 * C13451234 - m.m01 * C03451234 + m.m03 * C01451234 - m.m04 * C01351234 + m.m05 * C01341234);
	newM.m30 = det * -(m.m10 * C12452345 - m.m11 * C02452345 + m.m12 * C01452345 - m.m14 * C01252345 + m.m15 * C01242345);
	newM.m31 = det * (m.m00 * C12452345 - m.m01 * C02452345 + m.m02 * C01452345 - m.m04 * C01252345 + m.m05 * C01242345);
	newM.m32 = det * -(m.m00 * C12451345 - m.m01 * C02451345 + m.m02 * C01451345 - m.m04 * C01251345 + m.m05 * C01241345);
	newM.m33 = det * (m.m00 * C12451245 - m.m01 * C02451245 + m.m02 * C01451245 - m.m04 * C01251245 + m.m05 * C01241245);
	newM.m34 = det * -(m.m00 * C12451235 - m.m01 * C02451235 + m.m02 * C01451235 - m.m04 * C01251235 + m.m05 * C01241235);
	newM.m35 = det * (m.m00 * C12451234 - m.m01 * C02451234 + m.m02 * C01451234 - m.m04 * C01251234 + m.m05 * C01241234);
	newM.m40 = det * (m.m10 * C12352345 - m.m11 * C02352345 + m.m12 * C01352345 - m.m13 * C01252345 + m.m15 * C01232345);
	newM.m41 = det * -(m.m00 * C12352345 - m.m01 * C02352345 + m.m02 * C01352345 - m.m03 * C01252345 + m.m05 * C01232345);
	newM.m42 = det * (m.m00 * C12351345 - m.m01 * C02351345 + m.m02 * C01351345 - m.m03 * C01251345 + m.m05 * C01231345);
	newM.m43 = det * -(m.m00 * C12351245 - m.m01 * C02351245 + m.m02 * C01351245 - m.m03 * C01251245 + m.m05 * C01231245);
	newM.m44 = det * (m.m00 * C12351235 - m.m01 * C02351235 + m.m02 * C01351235 - m.m03 * C01251235 + m.m05 * C01231235);
	newM.m45 = det * -(m.m00 * C12351234 - m.m01 * C02351234 + m.m02 * C01351234 - m.m03 * C01251234 + m.m05 * C01231234);
	newM.m50 = det * -(m.m10 * C12342345 - m.m11 * C02342345 + m.m12 * C01342345 - m.m13 * C01242345 + m.m14 * C01232345);
	newM.m51 = det * (m.m00 * C12342345 - m.m01 * C02342345 + m.m02 * C01342345 - m.m03 * C01242345 + m.m04 * C01232345);
	newM.m52 = det * -(m.m00 * C12341345 - m.m01 * C02341345 + m.m02 * C01341345 - m.m03 * C01241345 + m.m04 * C01231345);
	newM.m53 = det * (m.m00 * C12341245 - m.m01 * C02341245 + m.m02 * C01341245 - m.m03 * C01241245 + m.m04 * C01231245);
	newM.m54 = det * -(m.m00 * C12341235 - m.m01 * C02341235 + m.m02 * C01341235 - m.m03 * C01241235 + m.m04 * C01231235);
	newM.m55 = det * (m.m00 * C12341234 - m.m01 * C02341234 + m.m02 * C01341234 - m.m03 * C01241234 + m.m04 * C01231234);
	return newM;
}
