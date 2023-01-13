#ifndef OMEGA_H_ESMFWRAPPER_H
#define OMEGA_H_ESMFWRAPPER_H
extern "C" void esmfInit();
extern "C" void esmfFinalize();
extern "C" void esmfTestMesh();
extern "C" void esmfLoadMesh(char* meshFileName, int strlen);
extern "C" void esmfGetMeshVtxIds(int* coords); //TODO use GO
extern "C" void esmfGetMeshVtxCoords(double* coords);
extern "C" void esmfGetMeshElemVerts(int* elemVerts);
#endif
