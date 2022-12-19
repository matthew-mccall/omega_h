#ifndef OMEGA_H_ESMFWRAPPER_H
#define OMEGA_H_ESMFWRAPPER_H
extern "C" void esmfInit();
extern "C" void esmfFinalize();
extern "C" void esmfTestMesh();
extern "C" void esmfGetMeshVtxCoords(double* coords);
#endif
