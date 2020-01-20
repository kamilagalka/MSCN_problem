#pragma once
#include <iostream>
#include "ErrorCodes.h"
template<typename  T>class Matrix
{
public:
	Matrix();
	Matrix(int sizeX, int sizeY);
	~Matrix();

	void changeSizeX(int newSizeX, int& errCode);
	void changeSizeY(int newSizeY, int& errCode);

	void setAt(int x, int y, T newValue, int& errCode);
	int getSizeX() { return sizeX; }
	int getSizeY() { return sizeY; }

	T getValueAt(int x, int y, int& errCode);

	void printMatrix();

private:
	int sizeX;
	int sizeY;
	T **table;

};

template<typename T>
inline Matrix<T>::Matrix()
{
	table = new T*[1];
	table[0] = new T[1];
	sizeX = 1;
	sizeY = 1;
}

template<typename T>
inline Matrix<T>::Matrix(int newSizeX, int newSizeY)
{

	if (newSizeX > 0 && newSizeY > 0) {
		sizeX = newSizeX;
		sizeY = newSizeY;

		table = new T*[sizeX];
		for (int i = 0; i < sizeX; i++) {
			table[i] = new T[sizeY];
		}

		for (int i = 0; i < sizeX; i++) {
			for (int j = 0; j < sizeY; j++) {
				table[i][j] = NULL;
			}
		}
	}
}

template<typename T>
inline Matrix<T>::~Matrix()
{
	for (int i = 0; i < sizeX; i++) {
		delete[] table[i];
	}
	delete[] table;
}

template<typename T>
inline void Matrix<T>::changeSizeX(int newSizeX, int& errCode)
{
	if (newSizeX <= 0) {
		errCode = NON_POSITIVE_VAL;
	}
	else {
		int smallerSize = sizeX < newSizeX ? sizeX : newSizeX;

		T **tempTab = new T*[newSizeX];

		for (int i = 0; i < smallerSize; i++) {
			tempTab[i] = table[i];
		}
		for (int i = smallerSize; i < newSizeX; i++) {
			tempTab[i] = new T[sizeY];

		}

		delete table;

		table = tempTab;

		sizeX = newSizeX;
	}


}

template<typename T>
inline void Matrix<T>::changeSizeY(int newSizeY, int& errCode)
{
	if (newSizeY <= 0) {
		errCode = NON_POSITIVE_VAL;
	}

	else {
		int smallerSize = sizeY < newSizeY ? sizeY : newSizeY;

		T **tempTab = new T*[sizeX];
		for (int i = 0; i < sizeX; i++) {
			tempTab[i] = new T[newSizeY];
		}

		for (int i = 0; i < sizeX; i++) {
			for (int j = 0; j < smallerSize; j++) {
				tempTab[i][j] = table[i][j];
			}
		}

		for (int i = 0; i < sizeX; i++) {
			delete table[i];
		}

		delete table;

		table = tempTab;

		sizeY = newSizeY;
	}


}

template<typename T>
inline void Matrix<T>::setAt(int x, int y, T value, int & errCode)
{
	if (x < 0 || x >= sizeX || y < 0 || y >= sizeY) {
		errCode = INVALID_OFFSET;
	}
	else {
		table[x][y] = value;
	}
}

template<typename T>
inline T Matrix<T>::getValueAt(int x, int y, int& errCode)
{
	if (x < 0 || x >= sizeX || y < 0 || y >= sizeY) {
		errCode = INVALID_OFFSET;
		return NULL;
	}
	return table[x][y];
}

template<typename T>
inline void Matrix<T>::printMatrix()
{
	for (int i = 0; i < sizeX; i++) {
		for (int j = 0; j < sizeY; j++) {
			std::cout << table[i][j] << "  ";
		}
		std::cout << std::endl;
	}
}
