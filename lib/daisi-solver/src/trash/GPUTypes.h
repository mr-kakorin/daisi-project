#ifndef GPUOBJECTS_H
#define GPUOBJECTS_H
#include "GPUTypes.h"
#include <cuda_runtime.h>
#include <iostream>
#include <vector>
namespace GPUTypes
{
class string
{
    char* data;
    int   DataSize;

  public:
    __device__ __host__ char& operator[](int i)
    {
        return data[i];
    }
    string(){};
    string(std::string str)
    {
        DataSize = str.size();
        cudaMalloc((void**)&data, DataSize * sizeof(char));
        cudaMemcpy(data, &str[0], DataSize * sizeof(char), cudaMemcpyHostToDevice);
    };
    string* AllocateOnGPU()
    {
        string* result;
        cudaMalloc((void**)&result, sizeof(string));
        cudaMemcpy(result, this, sizeof(string), cudaMemcpyHostToDevice);
        return result;
    };
    __device__ __host__ int size()
    {
        return DataSize;
    };
};
template <class T> class vector
{
    T*  data;
    int DataSize;

  public:
    __device__ __host__ T& operator[](int i)
    {
        return data[i];
    }
    vector(){};
    template <class TT> vector(std::vector<TT>& obj)
    {
        DataSize = obj.size();

        T* dataTmp = new T[DataSize];

        for (int i     = 0; i < DataSize; i++)
            dataTmp[i] = obj[i];

        cudaMalloc((void**)&data, DataSize * sizeof(T));
        cudaMemcpy(data, &dataTmp[0], DataSize * sizeof(T), cudaMemcpyHostToDevice);
        delete[] dataTmp;
    };

    vector<T>* AllocateOnGPU()
    {
        vector<T>* result;
        cudaMalloc((void**)&result, sizeof(vector<T>));
        cudaMemcpy(result, this, sizeof(vector<T>), cudaMemcpyHostToDevice);
        return result;
    };
    __device__ __host__ int size()
    {
        return DataSize;
    };
};

class Dmath::imat;
class imat
{
  private:
    vector<int> data;
    int         nrow;
    int         ncol;
    int         size;

  public:
    imat(Dmath::imat& obj)
    {
        data = obj.GetData();
        ncol = obj.GetNcol();
        nrow = obj.GetNrow();
        size = data.size();
    }
    /*	int getElem(unsigned int base, unsigned int d1, unsigned int d2)
            {
                    return data[base + d1 + d2*nrow];
            };
            int GetNrow()
            {
                    return nrow;
            };
            int ArSize()
            {
                    return int(data.size());
            };
            imat operator * (int right)
            {
                    imat result = *this;
                    for (int i = 0; i < size; i++)
                            result.data[i] = result.data[i] * right;
                    return result;
            };
            int& imat::operator() (int i, int j)
            {
                    return data[i + j*nrow];
            };
            imat & operator= (const imat & right)
            {
                    if (this == &right) {
                            return *this;
                    }
                    data = right.data;
                    nrow = right.nrow;
                    ncol = right.ncol;
                    size = right.size;
                    return *this;
            };*/
};
};
#endif