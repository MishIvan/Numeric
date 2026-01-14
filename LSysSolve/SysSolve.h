#pragma once
#include <iostream>
#include <chrono>
#include <cstdlib>
#include <string>
#include <Windows.h>
#include <WinUser.h>

#define MAX_ITERATION_NUMBER 60000

double rand_range(double min, double max);
double norm(const double* v, int n);
bool IsSymmetric(const double* A, int n);
void TransformLinearSystem(const double* A, const double* b, double* Anorm, double* bet, int n);

double Gauss(const double* a, const double* b, double* x, int size);
void CompactSchemeSolve(const double* A, const double* b, double* x, int n);
void QRDecompositionSolve(const double* A, const double* b, double* x, int n);
void LUDecompositionSolve(const double* A, const double* b, double* x, int n);
void LLTDecompositionSolve(const double* A, const double* b, double* x, int n);
bool Relaxation(const double* A, const double* b, double* x, int n, double omega);
void RotationSolve(const double* Asrc, const double* b, double* x, int n);