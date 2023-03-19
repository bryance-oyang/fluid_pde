#ifndef ARRAY_H
#define ARRAY_H

#include <cstring>
#include <utility>

#define MAX_ARRAY_RANK 5

/** contiguous in memory for cache efficiency accessed with A(i,j,k) */
template<typename T> class Array {
public:
	int rank;
	int n[MAX_ARRAY_RANK];
	int len;
	T *data = nullptr;

	void alloc()
	{
		data = new T[len];
	}
	void free()
	{
		if (data) {
			delete[] data;
			data = nullptr;
		}
	}
	size_t bytes() const
	{
		return len * sizeof(T);
	}
	void attach_reference(const Array &other)
	{
		len = other.len;
		rank = other.rank;
		memcpy(n, other.n, MAX_ARRAY_RANK * sizeof(n[0]));

		data = other.data;
	}
	void detach_reference()
	{
		len = 0;
		rank = 0;
		for (int i = 0; i < MAX_ARRAY_RANK; i++) {
			n[i] = 0;
		}
		data = nullptr;
	}
	void fill(const T &value)
	{
		for (int i = 0; i < len; i++) {
			data[i] = value;
		}
	}

	Array()
	: rank{0},
	n{0, 0, 0, 0, 0},
	len{0} {}

	Array(int n0)
	: rank{1},
	n{n0, 1, 1, 1, 1},
	len{n0} { alloc(); }

	Array(int n0, int n1)
	: rank{2},
	n{n0, n1, 1, 1, 1},
	len{n0*n1} { alloc(); }

	Array(int n0, int n1, int n2)
	: rank{3},
	n{n0, n1, n2, 1, 1},
	len{n0*n1*n2} { alloc(); }

	Array(int n0, int n1, int n2, int n3)
	: rank{4},
	n{n0, n1, n2, n3, 1},
	len{n0*n1*n2*n3} { alloc(); }

	Array(int n0, int n1, int n2, int n3, int n4)
	: rank{5},
	n{n0, n1, n2, n3, n4},
	len{n0*n1*n2*n3*n4} { alloc(); }

	~Array() { free(); }

	friend void swap(Array &first, Array &second)
	{
		using std::swap;
		swap(first.len, second.len);
		swap(first.rank, second.rank);
		swap(first.n, second.n);
		swap(first.data, second.data);
	}

	Array(const Array &other)
	{
		len = other.len;
		rank = other.rank;
		memcpy(n, other.n, MAX_ARRAY_RANK * sizeof(n[0]));

		free();
		alloc();
		memcpy(data, other.data, other.bytes());
	}

	Array(Array &&other)
	{
		swap(*this, other);
	}

	Array &operator=(Array other)
	{
		swap(*this, other);
		return *this;
	}

	void copy_data_from(Array &other)
	{
		memcpy(data, other.data, other.bytes());
	}

	T &operator()(int i0)
	{ return data[i0]; }
	T operator()(int i0) const
	{ return data[i0]; }

	T &operator()(int i0, int i1)
	{ return data[i0*n[1] + i1]; }
	T operator()(int i0, int i1) const
	{ return data[i0*n[1] + i1]; }

	T &operator()(int i0, int i1, int i2)
	{ return data[(i0*n[1] + i1)*n[2] + i2]; }
	T operator()(int i0, int i1, int i2) const
	{ return data[(i0*n[1] + i1)*n[2] + i2]; }

	T &operator()(int i0, int i1, int i2, int i3)
	{ return data[((i0*n[1] + i1)*n[2] + i2)*n[3] + i3]; }
	T operator()(int i0, int i1, int i2, int i3) const
	{ return data[((i0*n[1] + i1)*n[2] + i2)*n[3] + i3]; }

	T &operator()(int i0, int i1, int i2, int i3, int i4)
	{ return data[(((i0*n[1] + i1)*n[2] + i2)*n[3] + i3)*n[4] + i4]; }
	T operator()(int i0, int i1, int i2, int i3, int i4) const
	{ return data[(((i0*n[1] + i1)*n[2] + i2)*n[3] + i3)*n[4] + i4]; }
};

#endif /* ARRAY_H */
