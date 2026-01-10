#pragma once
#include <iostream>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <utility>

double rand_range(double min, double max);
double norm(const double* v, int n);
void TransformLinearSystem(const double* A, const double* b, double* Anorm, double* bet, int n);

double Gauss(const double* a, const double* b, double* x, int size);
void CompactSchemeSolve(const double* A, const double* b, double* x, int n);
void QRDecompositionSolve(const double* A, const double* b, double* x, int n);
void LUDecompositionSolve(const double* A, const double* b, double* x, int n);
void LLTDecompositionSolve(const double* A, const double* b, double* x, int n);