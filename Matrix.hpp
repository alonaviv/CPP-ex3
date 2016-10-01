//Matrix.hpp

#warning TODO unload long functions to non inline


# include <vector>
# include <stdexcept>
# include <iostream>
# include <algorithm>


const int ZERO_ARGUMENT = 0;
const int DEFAULT_ROW_SIZE = 1;
const int DEFAULT_COLUMN_SIZE = 1;

const int ITERATOR_STARTING_ROW = 0;
const int ITERATOR_STARTING_COL = 0;

const char CELL_SEPERATOR = '\t';

/** vector addition overloader. returns a new vector in which
 * each cell contains the addition of the two corresponding t 
 * objects in each vector. 
 * throws an excpetion if both vectors are not of the same size
 */
template<class t>
std::vector<t> operator+(const std::vector<t>& a, const std::vector<t>& b)
{
	if(a.size() != b.size())
	{
		throw std::length_error("trying to add to vectors of different lengths");
	}

	std::vector<t> result;
	for(unsigned int i = 0; i < a.size(); i++)
	{
		result.push_back(a[i] + b[i]);
	}
	return std::move(result);
}


/** vector substraction overloader. returns a new vector in which
 * each cell contains the substraction of the two corresponding
 * objects in each vector. 
 * throws an excpetion if both vectors are not of the same size
 */
template<class t>
std::vector<t> operator-(const std::vector<t>& a, const std::vector<t>& b)
{
	if(a.size() != b.size())
	{
		throw std::length_error("trying to add to vectors of different lengths");
	}

	std::vector<t> result;
	for(unsigned int i = 0; i < a.size(); i++)
	{
		result.push_back(a[i] - b[i]);
	}
	return std::move(result);
}


/**
 * This template class implements a container for a mathematical
 * matrix structure, complete with several types of constructors, various
 * mathematical operations, and returns an iterator for the container.
 * Implemented internally using a two dimensional vector
 *
 * IMPORTANT NOTE: Since Matrix only holds single data members, with
 * no memory management needed, it doesn't make use of the rule of five - 
 * it uses the default dtor, copy ctor, move assignment and copy assignment
 */
template<class T>
class Matrix
{
private:
	unsigned int _numberOfRows;
	unsigned int _numberOfCols;
	std::vector<std::vector<T>> _cells;

	/**
	 * Returns true iff given matrix has the same dimension as this
	 */
	bool _validateSameDimensions(const Matrix& other)
	{
		return _numberOfRows == other._numberOfRows and _numberOfCols == other._numberOfCols;
	}

	/**
	 * Receives two matrices, a row index for the first, and
	 * a column index for the second. Returns the dot product
	 * of the two vectors (in the linear algebra sense), of the row
	 * and the column. Throws an exception if the a's row size isn't equal
	 * to b's column size.
	 */
	static T s_matrixDotProduct(const Matrix& a, unsigned int aRowIndex, const Matrix& b,
			unsigned int bColIndex)
	{
		if(a._numberOfCols != b._numberOfRows)

		{
			std::cout << "a.number of cols is " << a._numberOfCols << std::endl<<"b.number of rows" << b._numberOfRows << std::endl << " aRowIndex" << aRowIndex << std::endl << "bColIndex" <<bColIndex << std::endl;
			throw std::length_error("Error computing dot product - incompatible sizes");
		}

		T sum(ZERO_ARGUMENT);
		for(unsigned int i = 0; i < a._numberOfCols; i++)
		{
			sum = sum + (a._cells[aRowIndex][i] * b._cells[i][bColIndex]);
		}
		return std::move(sum);
	}
	

public:
	
	/** Default ctor, creates a 1X1 matrix containing 
	 * T's zero element (created by calling C's ctor with 0
	 * as a parameter */
	Matrix(): Matrix(DEFAULT_ROW_SIZE, DEFAULT_COLUMN_SIZE){}


	/** Ctor that receives the dimensions of the matrix and creates a 
	 * corresponding matrix with the zero element of T in each cell.
	 */
	Matrix(unsigned int rows, unsigned int cols): _numberOfRows(rows), _numberOfCols(cols),
		   _cells(_numberOfRows, std::vector<T>(_numberOfCols, T(ZERO_ARGUMENT))){}
	 

	/** Ctor that receives the dimensions of a matrix, and a vector of data
	 * of type T. Create a matrix that holds the data in the vector, with
	 * each vector item being a single cell. Order is the same as matrix
	 * iteration (see constIterator object)
	 * Note: This ctor doesn't call Matrix(int,int), so not to create the
	 * vector already loaded with zero T objects. Instead, objects are
	 * inserted into the vector one by one.
	 * Throws exception if the vector's size isn't rows*cols.
	 */
	Matrix(unsigned int rows, unsigned int cols, const std::vector<T>& cells):
		_numberOfRows(rows), _numberOfCols(cols)
	{
		if(cells.size() != rows * cols)
		{
			throw std::out_of_range("Given vector size doesn't match matrix's dimensions");
		}

		for(unsigned int i = 0; i < _numberOfRows; i++)
		{
			 typename std::vector<T>::const_iterator rowStart = cells.cbegin() + i*_numberOfCols ;
			_cells.push_back(std::vector<T>(rowStart, rowStart + _numberOfCols));
		}
	}


	/** Ctor that receives the dimensions of a matrix and a vector of
	 * vectors (the exact data structure used by the matrix)
	 */
	Matrix(unsigned int rows, unsigned int cols, std::vector<std::vector<T>> cells):
		_numberOfRows(rows), _numberOfCols(cols), _cells(cells){}


	/* Performs the addition operation on two matrices of the same size. 
	 * Throws an exception if the two matrices are not of the same size
	 */
	Matrix operator+(const Matrix<T>& other) const
	{
		
		if(not _validateSameDimensions(other))
		{
			throw std::length_error("Trying to add matrices of different dimensions");
		}

		std::vector<std::vector<T>> newCells;
		for(unsigned int i = 0; i < _numberOfRows; i++)
		{
			newCells.push_back(_cells[i] + other._cells[i]);
		}
		return std::move(Matrix(_numberOfRows, _numberOfCols, std::move(newCells)));
	}


	/* Performs the subtraction operation on two matrices of the same size. 
	 * Throws an exception if the two matrices are not of the same size
	 */
	Matrix operator-(const Matrix<T>& other) const
	{
		
		if(not _validateSameDimensions(other))
		{
			throw std::length_error("Trying to add matrices of different dimensions");
		}

		std::vector<std::vector<T>> newCells;
		for(unsigned int i = 0; i < _numberOfRows; i++)
		{
			newCells.push_back(_cells[i] - other._cells[i]);
		}
		return std::move(Matrix(_numberOfRows, _numberOfCols, std::move(newCells)));
	}


	/* Performs the multiplication operation on two matrices, according to
	 * the linear algebra definition of multiplication.
	 * Throws an exception if the two matrices are not of compatible
	 * multiplication dimensions.
	 */
	Matrix operator*(const Matrix<T>& other) const 
	{
		if(_numberOfCols != other._numberOfRows)
		{
			throw std::length_error("Trying to add matrices of different dimensions");
		}

		std::vector<std::vector<T>> cells;
		
		std::vector<T> row;
		for(unsigned int i = 0; i < _numberOfRows; i++)
		{
			for(unsigned int j = 0; j < other._numberOfCols; j++)
			{
				row.push_back(std::move(s_matrixDotProduct(*this, i, other, j)));
			}
			cells.push_back(move(std::move(row)));
		}
		return std::move(Matrix(_numberOfRows, other._numberOfCols, std::move(cells)));
	}


	/* Performs a transposition of the matrix. Returns a new, transformed 
	 * copy, without modifying the original.*/
	 
	Matrix trans() const 
	{
		std::vector<std::vector<T>> cells;
		
		std::vector<T> row;
		for(unsigned int i = 0; i < _numberOfCols; i++)
		{
			for(unsigned int j = 0; j < _numberOfRows; j++)
			{
				row.push_back(_cells[j][i]);
			}
			cells.push_back(move(std::move(row)));
		}
		return std::move(Matrix(_numberOfCols, _numberOfRows, std::move(cells)));
	}


	/**
	 * Equality operator. Returns true iff all items within the matrix are
	 * equal
	 */
	bool operator==(const Matrix<T>& other) const 
	{
		return _cells == other._cells;
	}


	/**
	 * Inequality operator. Returns false iff all items within the matrix are
	 * equal
	 */
	bool operator!=(const Matrix<T>& other) const 
	{
		return not (*this == other);
	}


	/** 
	 * returns a const ref of the contents of a cell according to given indexes.
	 * throws an exception if the given indexes exceeds the matrix's boundaries.
	 */
	const T& operator()(unsigned int rowIndex, unsigned int colIndex) const
	{
		if(rowIndex >= _numberOfRows or colIndex >= _numberOfCols)
		{
			throw std::out_of_range("requested index exceeds matrix's range");
		}
		return _cells[rowIndex][colIndex];
	}

	
	/** 
	 * returns a non-const ref of the contents of a cell according to given indexes.
	 * throws an exception if the given indexes exceeds the matrix's boundaries.
	 */
# warning TODO Check that both constant/non-constant are called when T is a object
	T& operator()(unsigned int rowIndex, unsigned int colIndex) 
	{
		if(rowIndex >= _numberOfRows or colIndex >= _numberOfCols)
		{
			throw std::out_of_range("requested index exceeds matrix's range");
		}
		return _cells[rowIndex][colIndex];
	}


	/** Returns true iff the matrix is a square */
	bool isSquareMatrix() const
	{
		return _numberOfRows == _numberOfCols;
	}


	/** 
	 * Returns the number of rows in the matrix
	 */
	unsigned int rows() const
	{
		return _numberOfRows;
	}


	/** 
	 * Returns the number of columns in the matrix.
	 */
	unsigned int cols() const
	{
		return _numberOfCols;
	}

	
	/** Outputs a matrix to a given ostream. Each value of a 
	 * column is separated by a \t, and each row is
	 * separated by a new line */
	friend std::ostream& operator<<(std::ostream &os, const Matrix& matrix)
	{
		for(auto i = matrix._cells.begin(); i != matrix._cells.end(); i++)
		{
			for(auto j = i -> begin(); j != i -> end(); j++)
			{
				os << *j << CELL_SEPERATOR;
			}
			os << std::endl;
		}
		return os;
	}	

#warning TODO add -- oper plus any other standard iter oper I need
	/** Bidirectional constant iterator object. Iterates over the matrix's values
	 * by traversing the matrix row by row, starting from the left.
	 * The end of the matrix is defined as being in location
	 * [rowSize][0] - meaning, one cell after the last cell in the
	 * matrix.
	 */
	class constIterator
	{
	private:
		unsigned int _rowLocation;
		unsigned int _colLocation;
		const Matrix& _matrix;

	public:
		/**
		 * Constructor. Receives a pointer to the starting 
		 * location of the iterator withing the matrix.
		 */
		constIterator(const unsigned int startingRowLocation,
				const unsigned int startingColLocation, 
				const Matrix& matrix): _rowLocation(startingRowLocation),
		_colLocation(startingColLocation), _matrix(matrix){}
				
		/**
		 * Dereferencing operator, returns a refernece the T object that the
		 * iterator is pointing to
		 */
		const T& operator*() const
		{
			return _matrix._cells[_rowLocation][_colLocation];
		}

		/**
		 * Arrow operator. Returns the raw pointer to the current 
		 * T object
		 */
		const T* operator->() const
		{
			return &_matrix.cells[_rowLocation][_colLocation];
		}

		/**
		 * prefix increment iterator.
		 * throws an exception if iterator is already at the last
		 * location.
		 */
		constIterator& operator++()
		{
			if (_rowLocation ==  _matrix._numberOfRows and _colLocation == ITERATOR_STARTING_COL)
			{
				throw std::out_of_range("iterator is already at matrix's end");
			}

			if(_colLocation != _matrix._numberOfCols - 1)
			{
				_colLocation++;
			}
			else
			{
				_colLocation = ITERATOR_STARTING_COL;
				_rowLocation ++;
			}
			return *this;
		}
		

		/**
		 * prefix decrement iterator.
		 * throws an exception if iterator is already at the first 
		 * location.
		 */
		constIterator& operator--()
		{
			if (_rowLocation == ITERATOR_STARTING_ROW and _colLocation == ITERATOR_STARTING_COL)
			{
				throw std::out_of_range("iterator is already at matrix's start");
			}

			if(_colLocation != ITERATOR_STARTING_COL)
			{
				_colLocation--;
			}
			else
			{
				_colLocation = _matrix._numberOfCols - 1;
				_rowLocation --;
			}
			return *this;
		}


		/**
		 * Postfix Increment iterator.
		 * Throws an exception if iterator is already at the last
		 * location.
		 */
		constIterator operator++(int)
		{
			constIterator currentIterator(*this);
			++(*this);
			return currentIterator;
		}

		/**
		 * Postfix decrement iterator.
		 * Throws an exception if iterator is already at the last
		 * location.
		 */
		constIterator operator--(int)
		{
			constIterator currentIterator(*this);
			--(*this);
			return currentIterator;
		}

		/**
		 * Equality Comparison between two iterators, by their location 
		 * in the matrix
		 */
		bool operator==(const constIterator& other) const
		{
			return _rowLocation == other._rowLocation
				and _colLocation == other._colLocation;
		}

		/**
		 * Unequality Comparison between two iterators, by their location
		 * in the matrix 
		 */
		bool operator!=(const constIterator& other) const
		{
			return not (*this == other);
		}
	};

	/** Returns an iterator to the start of the matrix */
	constIterator begin() const
	{
		return constIterator(ITERATOR_STARTING_ROW, ITERATOR_STARTING_COL, *this);
	}

	/** Returns an iterator to the end of the matrix. Defined 
	 * in the context as the cell after the last, i.e [rowSize][0]
	 */
	constIterator end() const
	{
		return constIterator(_numberOfRows, ITERATOR_STARTING_COL, *this);
	}

};


