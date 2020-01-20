#pragma once
#include "ErrorCodes.h"
#include <iostream>

template <typename T>class Array {

public:
	Array() { table = new T[1]; size = 1; };
	Array(int newSize) { table = new T[newSize]; size = newSize; }
	~Array() { delete table; }

	void changeSize(int newSize, int& errCode);
	void setAt(int index, T value, int& errCode);
	int getSize() { return size; }
	T getValueAt(int index, int& errCode);

	T* getTable() { return table; }

	void printArray();

private:
	int size;
	T *table;
};



template<typename T>
inline void Array<T>::changeSize(int newSize, int & errCode)
{
	if (newSize <= 0) {
		errCode = NON_POSITIVE_VAL;
	}
	
	else {
		T *tempTab = new T[newSize];

		int smallerSize = size < newSize ? size : newSize;

		for (int i = 0; i < smallerSize; i++) {
			tempTab[i] = table[i];
		}
		delete table;
		table = tempTab;
		size = newSize;
	}
}

template<typename T>
inline void Array<T>::setAt(int index, T value, int & errCode)
{
	if (index < 0 || index >= size) {
		errCode = INVALID_OFFSET;
	}
	else {
		table[index] = value;
	}

}

template<typename T>
inline T Array<T>::getValueAt(int index, int & errCode)
{
	if (index < 0 || index >= size) {
		errCode = INVALID_OFFSET;
		return NULL;
	}
	return table[index];
}

template<typename T>
inline void Array<T>::printArray()
{
	for (int i = 0; i < size; i++) {
		std::cout << table[i] << " ";
	}
}
