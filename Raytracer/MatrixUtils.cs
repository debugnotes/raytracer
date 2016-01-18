using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace MatrixUtils
{

	public class Vector
	{
		protected double[] arrVector = null;

		/// <summary>
		/// Constructor.
		/// </summary>
		public Vector(int iSize)
		{
			Resize(iSize, true);
		}

		/// <summary>
		/// Constructor for a vector of size 3, or a 3D vector.
		/// </summary>
		/// <param name="dX"></param>
		/// <param name="dY"></param>
		/// <param name="dZ"></param>
		public Vector(double dX, double dY, double dZ)
		{
			Resize(3, true);
			arrVector[0] = dX;
			arrVector[1] = dY;
			arrVector[2] = dZ;
		}

		/// <summary>
		/// Copy constructor.
		/// </summary>
		/// <param name="?"></param>
		public Vector(Vector vectSource)
		{
			if (vectSource == null)
				throw new Exception("Error in Vector.Vector() - source vector cannot be null.");
			Resize(vectSource.arrVector.Length, false);
			for (int iIndex = 0; iIndex < vectSource.arrVector.Length; iIndex++)
			{
				arrVector[iIndex] = vectSource.arrVector[iIndex];
			}
		}

		/// <summary>
		/// Creates a matrix with the specified dimensions. If 
		/// bFastCreate, then the matrix is not initialized, and
		/// the contents of the matrix are undefined.
		/// </summary>
		public Vector(int iSize, bool bFastCreate)
		{
			if (bFastCreate)
				arrVector = Create1DArray(iSize);
			else
				Resize(iSize, true);
		}

		/// <summary>
		/// Gets or sets the value of a matrix entry at [iIndex].
		/// </summary>
		public double this[int iIndex]
		{
			get
			{
				return arrVector[iIndex];
			}
			set
			{
				arrVector[iIndex] = value;
			}
		}

		/// <summary>
		/// Returns true if the matrix is null, false otherwise.
		/// </summary>
		public bool IsNull
		{
			get
			{
				if (arrVector == null)
					return true;
				return false;
			}
		}
		/// <summary>
		/// Gets the unit length of the vector. Note that this is
		/// NOT the same as Size, which returns the number of
		/// elements in the Vector.
		/// </summary>
		public double Length
		{
			get
			{
				if (arrVector == null)
					return 0;
				int iSize = Size;
				double dAccumulator = 0;
				for (int iIndex = 0; iIndex < iSize; iIndex++)
				{
					dAccumulator += (arrVector[iIndex] * arrVector[iIndex]);
				}
				return Math.Sqrt(dAccumulator);
			}
		}

		/// <summary>
		/// Gets the size of the vector (i.e. number of elements)
		/// </summary>
		public int Size
		{
			get
			{
				if (arrVector == null)
					return 0;
				return arrVector.Length;
			}
			set
			{
				Resize(value, false);
			}
		}

		public static Vector operator+ (Vector vect1, Vector vect2)
		{
			if (vect1 == null || vect2 == null)
				throw new NullReferenceException();
			if (vect1.arrVector == null || vect2.arrVector == null)
				throw new Exception("Addition of null vectors is not permitted.");
			if (vect1.arrVector.Length != vect2.arrVector.Length)
				throw new Exception("Cannot add - vectors are not of the same size.");
			Vector vectResult = new Vector(vect1.arrVector.Length, true);
			for (int iIndex = 0; iIndex < vect1.arrVector.Length; iIndex++)
			{
				vectResult.arrVector[iIndex] = vect1.arrVector[iIndex] + vect2.arrVector[iIndex];
			}
			return vectResult;
		}


		public static Vector operator- (Vector vect1, Vector vect2)
		{
			if (vect1 == null || vect2 == null)
				throw new NullReferenceException();
			if (vect1.arrVector == null || vect2.arrVector == null)
				throw new Exception("Subtraction of null vectors is not permitted.");
			if (vect1.arrVector.Length != vect2.arrVector.Length)
				throw new Exception("Cannot subtract - vectors are not of the same size.");
			Vector vectResult = new Vector(vect1.arrVector.Length, true);
			for (int iIndex = 0; iIndex < vect1.arrVector.Length; iIndex++)
			{
				vectResult.arrVector[iIndex] = vect1.arrVector[iIndex] - vect2.arrVector[iIndex];
			}
			return vectResult;
		}

		public static Vector operator* (double dVal, Vector vect)
		{
			return vect * dVal;
		}

		public static Vector operator* (Vector vect, double dVal)
		{
			if (vect == null)
				throw new NullReferenceException();
			if (vect.arrVector == null)
				throw new Exception("Scaling a null vector is not permitted.");
			Vector vectResult = new Vector(vect.arrVector.Length, true);
			for (int iIndex = 0; iIndex < vect.arrVector.Length; iIndex++)
			{
				vectResult.arrVector[iIndex] = dVal * vect.arrVector[iIndex];
			}
			return vectResult;
		}

		public static Vector operator* (Vector vect1, Vector vect2) {
			if (vect1 == null || vect2 == null)
				throw new NullReferenceException();
			if (vect1.arrVector == null || vect2.arrVector == null)
				throw new Exception("Multiplying by null vector is not permitted.");
			if (vect1.arrVector.Length != vect2.arrVector.Length)
				throw new Exception("Multiplying two vectors together by different lengths is not permitted.");
			Vector vectResult = new Vector(vect1.arrVector.Length, true);
			for (int iIndex = 0; iIndex < vect1.arrVector.Length; iIndex++)
			{
				vectResult.arrVector[iIndex] = vect1.arrVector[iIndex] * vect2.arrVector[iIndex];
			}
			return vectResult;
		}

		/// <summary>
		/// Calculates the angle betweeen two vectors. Angle returned is in degrees.
		/// </summary>
		public double Angle(Vector vect)
		{
			double dVal = DotProduct(vect) / (Length * vect.Length);
			return Math.Acos(dVal) * (180.0 / Math.PI);
		}

		/// <summary>
		/// Sets all elements of the vector zero.
		/// </summary>
		public void Clear()
		{
			if (arrVector == null)
				return;
			int iSize = Size;
			for (int iIndex = 0; iIndex < iSize; iIndex++)
			{
				arrVector[iIndex] = 0;
			}
		}

		/// <summary>
		/// Returns the dot product of this and the second provided vector,
		/// </summary>
		/// <param name="vect"></param>
		/// <returns></returns>
		public double DotProduct(Vector vect)
		{
			if (vect == null)
				return 0;
			if (vect.arrVector == null || arrVector == null)
				throw new Exception ("Error in Vector.DotProduct(). The vector cannot be a null vector.");
			if (vect.Size != Size)
				throw new Exception ("Error in Vector.DotProduct(). Both vectors cannot be of the same size.");
			double dResult = 0;
			for (int iIndex = 0; iIndex < arrVector.Length; iIndex++)
			{
				dResult += vect.arrVector[iIndex] * arrVector[iIndex]; 
			}
			return dResult;
		}

		/// <summary>
		/// Creates a new vector that is the unit vector for
		/// for this vector.
		/// </summary>
		public Vector Normalize()
		{
			if (arrVector == null)
				return new Vector(0);
			Vector vectNormalized = new Vector(arrVector.Length, true);
			double dUnitLength = Length;
			if (arrVector != null && dUnitLength > 0)
			{
				for (int iIndex = 0; iIndex < arrVector.Length; iIndex++)
				{
					vectNormalized[iIndex] = arrVector[iIndex] / dUnitLength;
				}
			}
			return vectNormalized;
		}



		/// <summary>
		/// Sets the size of the matrix to the specified value.
		/// All data that is in the current matrix will be copied to the 
		/// new matrix if possible.
		/// </summary>
		public void Resize(int iSize)
		{
			Resize(iSize, false);
		}

		/// <summary>
		/// Sets the size of the matrix to the specified value.
		/// All data that is in the current matrix will be copied to the 
		/// new matrix if possible, unless bClearData is set to true.
		/// </summary>
		public void Resize(int iSize, bool bClearData)
		{
			if (iSize < 1)
			{
				arrVector = null;
				return;
			}
			double[] arrNewVector = new double[iSize];
			int iOldSize = Size;
			for (int iIndex = 0; iIndex < iSize; iIndex++)
			{
				if (arrVector != null && bClearData == false &&
					iIndex < iOldSize)
				{
					arrNewVector[iIndex] = arrVector[iIndex];
				}
				else
				{
					arrNewVector[iIndex] = 0;
				}
			}
			// set the old vector to the new
			arrVector = arrNewVector;
		}

		/// <summary>
		/// Produces a string that describes the elements of the
		/// vector.
		/// </summary>
		public override String ToString()
		{
			if (arrVector == null)
				return null;
			StringBuilder sbOutput = new StringBuilder();
			int iSize = Size;
			sbOutput.AppendLine("Vector of size " + iSize);
			for (int iIndex = 0; iIndex < iSize; iIndex++)
			{
				sbOutput.Append(String.Format("{0:G4}", arrVector[iIndex]));
				if (iIndex < iSize - 1)
					sbOutput.Append(", ");
				sbOutput.AppendLine();
			}
			return sbOutput.ToString();
		}

		/// <summary>
		/// Creates a new vector that is the unit vector for
		/// for this vector. Exactly the same as Normalize().
		/// </summary>
		public Vector UnitVector()
		{
			return Normalize();
		}

		/// <summary>
		/// Adjusts the values of a vector so that they are within
		/// the limits imposed by Min and Max.
		/// </summary>
		/// <param name="dMin"></param>
		/// <param name="dMax"></param>
		/// <returns></returns>
		public Vector Limit(double dMin, double dMax)
		{
			if (arrVector == null)
				return new Vector(0);
			if (dMin > dMax)
			{
				double dTemp = dMin;
				dMin = dMax;
				dMax = dTemp;
			}
			Vector vectLimited = new Vector(arrVector.Length, true);
			for (int iIndex = 0; iIndex < vectLimited.Size; iIndex++)
			{
				if (arrVector[iIndex] < dMin)
					vectLimited[iIndex] = dMin;
				else if (arrVector[iIndex] > dMax)
					vectLimited[iIndex] = dMax;
				else
					vectLimited[iIndex] = arrVector[iIndex];
			}
			return vectLimited;
		}
		/// <summary>
		/// Sets all elements of the vector zero. Behaves exactly
		/// the same way as zero.
		/// </summary>
		public void Zero()
		{
			Clear();
		}

		/// <summary>
		/// Creates a 1D array of doubles, without the initialization
		/// overhead of Resize.
		/// </summary>
		protected double[] Create1DArray(int iSize)
		{
			if (iSize < 0)
				return null;
			return new double[iSize];
		}
	}; // class Vector

	/// <summary>
	/// A generic matrix class.
	/// </summary>
	public class Matrix
	{
		protected double[,] arrMatrix;
		int iRows, iColumns;

		/// <summary>
		/// Default constructor.
		/// </summary>
		public Matrix (int iRows, int iColumns)
		{
			Resize(iRows, iColumns, true);
		}

		/// <summary>
		/// Creates a matrix with the specified dimensions. If 
		/// bFastCreate, then the matrix is not initialized, and
		/// the contents of the matrix are undefined.
		/// </summary>
		public Matrix (int iRows, int iColumns, bool bFastCreate)
		{
			if (bFastCreate)
				arrMatrix = Create2DArray(iRows, iColumns);
			else
				Resize(iRows, iColumns, true);
		}

		/// <summary>
		/// Copy constructor.
		/// </summary>
		public Matrix (Matrix mat)
		{
			Assign(mat);
		}

		/// <summary>
		/// Copy constructor with a vector.
		/// </summary>
		public Matrix (Vector vect)
		{
			Assign(vect);
		}

		/// <summary>
		/// Gets the nummber of columns in the matrix, and sets
		/// the number of columns indirectly by calling Resize()
		/// </summary>
		public int Columns
		{
			get
			{
				if (arrMatrix == null)
					return 0;
				return iColumns;
			}
			set
			{
				Resize(Rows, value, false);
			}
		}

		/// <summary>
		/// Gets or sets the value of a matrix entry at [iRow, iCol]
		/// </summary>
		public double this[int iRow, int iCol]
		{
			get
			{
				return arrMatrix[iRow,iCol];
			}
			set
			{
				arrMatrix[iRow,iCol] = value;
			}
		}

		/// <summary>
		/// Returns true if the matrix is diagonal, which can only be true
		/// if the matrix is square and there exists only numbers down the main
		/// diagonal of the matrix.
		/// </summary>
		public bool IsDiagonal
		{
			get
			{
				if (iRows != iColumns)
					return false;
				bool bIsDiagonal = true;
				for (int iRow = 0; iRow < iRows && bIsDiagonal == true; iRow++)
				{
					for (int iCol = 0; iCol < iColumns; iCol++)
					{
						if (iCol == iRow)
							continue;
						if (arrMatrix[iRow, iCol] != 0)
						{
							bIsDiagonal = false;
							break;
						}
					}
				}				
				return false;
			}
		}

		/// <summary>
		/// Returns true if the matrix is the identity matrix, which must
		/// be a square matrix as well.
		/// </summary>
		public bool IsIdentity
		{
			get
			{
				if (arrMatrix == null)
					return false;
				if (iRows != iColumns)
					return false;
				bool bIsIdentity = true;
				for (int iRow = 0; iRow < iRows && bIsIdentity == true; iRow++)
				{
					for (int iCol = 0; iCol < iColumns; iCol++)
					{
						if (iCol == iRow && arrMatrix[iRow, iCol] != 1)
						{
							bIsIdentity = false;
							break;
						}
						else if (arrMatrix[iRow,iCol] != 0)
						{
							bIsIdentity = false;
							break;
						}
					}
				}				
				return bIsIdentity;
			}
		}

		/// <summary>
		/// Returns true if the matrix is invertible. The matrix
		/// MUST be square to be invertible.
		/// </summary>
		public bool IsInvertible
		{
			get
			{
				if (arrMatrix == null)
					return false;
				if (iRows == iColumns)
					return false;
				// 1. Do gauss jordan reduction on this matrix.
				Matrix matTemp = GaussJordonElimination();
				// 2. See if the resulting matrix on LHS is the identity
				matTemp.Resize(3, 3, false);
				if (matTemp.IsIdentity)
					return true;
				return false;
			}
		}

		/// <summary>
		/// Returns true if the matrix is null, false otherwise.
		/// </summary>
		public bool IsNull
		{
			get
			{
				if (arrMatrix == null)
					return true;
				return false;
			}
		}

		/// <summary>
		/// Returns true if the matrix is square (i.e. has same number
		/// of rows and columns), false otherwise.
		/// </summary>
		public bool IsSquare
		{
			get
			{
				if (Rows == Columns)
					return true;
				return false;
			}
		}

		/// <summary>
		/// Returns true if the matrix is a zero matrix (i.e. all
		/// entries are zero), false otherwise.
		/// </summary>
		public bool IsZero
		{
			get
			{
				if (arrMatrix == null)
					throw new Exception ("Error in Matrix.IsZero. " +
						"This matrix has a size of zero, so it is undefined.");
				bool bIsZero = true;	// assume true
				for (int iRow = 0; iRow < iRows && bIsZero; iRow++)
				{
					for (int iCol = 0; iCol < iColumns; iCol++)
					{
						if (arrMatrix[iRow, iCol] != 0)
						{
							bIsZero = false;
							break;
						}
					}
				}
				return bIsZero;
			}
		}

		/// <summary>
		/// Gets the nummber of rows in the matrix, and sets
		/// the number of rows indirectly by calling Resize()
		/// </summary>
		public int Rows
		{
			get
			{
				if (arrMatrix == null)
					return 0;
				return arrMatrix.Length;
			}
			set
			{
				Resize(value, Columns, true);
			}
		}

		/// <summary>
		/// Returns the trace of a matrix, which is the sum of all the
		/// diagonals in a square matrix. Throws an exception if the
		/// matrix is not square.
		/// </summary>
		public double Trace
		{
			get
			{
				if (arrMatrix == null)
					throw new Exception ("Error in Matrix.Trace. " +
						"This matrix has a size of zero, so trace is undefined.");
				int iRows = Rows;
				if (iRows != Columns)
					throw new Exception ("Error in Matrix.Trace. " +
						"This matrix is not square, so trace is undefined.");
				double dResult = 0;
				for (int iRow = 0; iRow < iRows; iRow++)
				{
					dResult += arrMatrix[iRow, iRow];
				}
				return dResult;
			}
		}

		/// <summary>
		/// Makes the matrix an exact duplicate of the provided matrix.
		/// </summary>
		public void Assign(Matrix mat)
		{
			int iRows = mat.iRows;
			int iColumns = mat.iColumns;
			arrMatrix = Create2DArray(iRows, iColumns);
			if (arrMatrix == null)
				return;
			for (int iRow = 0; iRow < iRows; iRow++)
			{
				for (int iCol = 0; iCol < iColumns; iCol++)
				{
					arrMatrix[iRow ,iCol] = mat.arrMatrix[iRow ,iCol];
				}
			}
		}

		/// <summary>
		/// Makes the matrix an exact duplicate of the provided vector.
		/// </summary>
		public void Assign(Vector vect)
		{
			int iColumns = vect.Size;
			arrMatrix = Create2DArray(1, iColumns);
			if (arrMatrix == null)
				return;
			for (int iCol = 0; iCol < iColumns; iCol++)
			{
				arrMatrix[0, iCol] = vect[iCol];
			}
		}

		/// <summary>
		/// Clears the matrix by setting all entries to zero.
		/// </summary>
		public void Clear()
		{
			if (arrMatrix == null)
				return;
			for (int iRow = 0; iRow < iRows; iRow++)
			{
				for (int iCol = 0; iCol < iColumns; iCol++)
				{
					arrMatrix[iRow, iCol] = 0;
				}
			}
		}

		/// <summary>
		/// Creates a copy of this matrix.
		/// </summary>
		public Matrix Clone()
		{
			return new Matrix(this);
		}

		/// <summary>
		/// Returns true if two matricies are equal, which means they
		/// both must be the same size, and both must have the same
		/// corresponding entries.
		/// </summary>
		public virtual bool Equals(Matrix mat)
		{
			if (mat == null)
				return false;
			if (arrMatrix == null && mat.arrMatrix == null)
				return true;
			if ((arrMatrix != null && mat.arrMatrix == null) ||
				(arrMatrix == null && mat.arrMatrix != null))
				return false;
			bool bIsEqual = false;		// assume it isn't
			if (mat.Rows == iRows && mat.Columns == iColumns)
			{
				bIsEqual = true;		// assume they are equal
				for (int iRow = 0; iRow < iRows && bIsEqual == true; iRow++)
				{
					for (int iCol = 0; iCol < iColumns; iCol++)
					{
						if (arrMatrix[iRow,iCol] != mat.arrMatrix[iRow,iCol])
						{
							bIsEqual = false;
							break;
						}
					}
				}
			}
			return bIsEqual;
		}

		/// <summary>
		/// Performs Gause Jordan Elimination on the matrix.
		/// </summary>
		public Matrix GaussJordonElimination()
		{
			if (arrMatrix == null)
				return new Matrix(0, 0);
			Matrix matResult = new Matrix(this);
			int iMajorColumn = 0;
			for (int iMajorRow = 0; iMajorRow < iRows; iMajorRow++)
			{
				bool bLeadingEntryFound = false;
				while (!bLeadingEntryFound && iMajorColumn < iColumns)
				{
					// step 1 - if the row has a leading zero in the
					// column with the same index, find a row
					// with a leading entry
					if (matResult.arrMatrix[iMajorRow, iMajorColumn] == 0)
					{
						double[] arrTemp = new double[iColumns];
						for (int iRow = iMajorRow + 1; iRow < iRows; iRow++)
						{
							if (matResult.arrMatrix[iRow, iMajorColumn] != 0)
							{
								// swap the current row with the major row
								for (int iCol = 0; iCol < iColumns; iCol++)
								{
									arrTemp[iCol] = matResult.arrMatrix[iMajorRow, iCol];
									matResult.arrMatrix[iMajorRow, iCol] = matResult.arrMatrix[iRow, iCol];
									matResult.arrMatrix[iRow, iCol] = arrTemp[iCol];
								}
								bLeadingEntryFound = true;
								break;
							}
						}
						// move on to the next column
						if (!bLeadingEntryFound)
							iMajorColumn++;
					}
					else
					{
						bLeadingEntryFound = true;
					}
				}
				// if we can't find a row with a non leading zero, then we're done
				if (bLeadingEntryFound == false)
					break;

				//Debug.WriteLine(ToString());
				// step 2 - multiply the entire current row by the inverse
				// of the leading element, so that the leading element has
				// a value of 1
				double dLeadingInverse = 1.0/matResult.arrMatrix[iMajorRow, iMajorColumn];
				for (int iCol = iMajorColumn; iCol < iColumns; iCol++)
				{
					matResult.arrMatrix[iMajorRow, iCol] *= dLeadingInverse;
				}
				//Debug.WriteLine(ToString());
				// step 3 - go through remaining rows, and add mutiples of
				// the current major row to make the leading entry in the
				// corresponding column 0 below the leading 1
				for (int iRow = iMajorRow + 1; iRow < iRows; iRow++)
				{
					if (matResult.arrMatrix[iRow, iMajorColumn] == 0)
						continue;
					double dZeroedValue = -1 * matResult.arrMatrix[iRow, iMajorColumn];
					// add the zeroed value as a multiple of the major row
					// to the current row
					for (int iCol = iMajorColumn; iCol < iColumns; iCol++)
					{
						matResult.arrMatrix[iRow, iCol] += 
							matResult.arrMatrix[iMajorRow, iCol]* dZeroedValue;
					}
				}
				//Debug.WriteLine(ToString());
				//Debug.WriteLine("Row" + iMajorRow + " complete.");
				iMajorColumn++;
			}
			// step 4 - now go from the bottom of the matrix upwards
			// to introduce zeros above the leading ones in each row
			// to produce a matrix in reduced row echelon format
			iMajorColumn = iColumns - 1;
			for (int iMajorRow = iRows - 1; iMajorRow >= 0; iMajorRow--)
			{
				// find the leading 1 in the row (reading from right to left)
				for (int iCol = 0; iCol <= iMajorColumn; iCol++)
				{
					if (matResult.arrMatrix[iMajorRow, iCol] != 0)
					{
						iMajorColumn = iCol;
						break;
					}
				}
				// if the major column was not found, then the row is
				// all zeros, so move on to the next row
				if (matResult.arrMatrix[iMajorRow, iMajorColumn] == 0)
					continue;
				for (int iRowAbove = iMajorRow-1; iRowAbove >= 0; iRowAbove--)
				{
					double dZeroedValue = -1 * matResult.arrMatrix[iRowAbove, iMajorColumn]; 
					for (int iCol = iColumns - 1; iCol >= iMajorColumn; iCol--)
					{
						matResult.arrMatrix[iRowAbove, iCol] += 
							matResult.arrMatrix[iMajorRow, iCol]* dZeroedValue;
					}
				}
			}
			return matResult;
		}

		/// <summary>
		/// Returns the inverse of a square matrix. Throws an exception
		/// if the matrix is null or not quare - use the IsIdentity
		/// property to determine if the matrix can be inverted. If the
		/// matrix is not invertible, a null matrix is returned. Use the
		/// IsNull property on the resultant matrix to determine if 
		/// the matrix has an inverse.
		/// </summary>
		public Matrix Inverse()
		{
			// check for a null matrix
			if (arrMatrix == null)
			{
				throw new Exception("Error in Matrix.Inverse(): " +
					"The matrix is null, and cannot be inverted.");
			}
			if (iRows != iColumns)
			{
				throw new Exception("Error in Matrix.Inverse(): " +
					"The matrix is not a square matrix, and cannot be inverted.");
			}
			Matrix matResult = new Matrix(this);
			// ok, resize the matrix to twice the number of columns
			//
			// A B C | 0 0 0
			// D E F | 0 0 0
			// G H I | 0 0 0
			//
			matResult.Resize(iRows, iColumns*2);
			// setup a new identity matrix on the RHS, e.g.:
			//
			// A B C | 1 0 0
			// D E F | 0 1 0
			// G H I | 0 0 1
			//
			for (int iRow = 0; iRow < iRows; iRow++)
			{
				matResult.arrMatrix[iRow, iColumns+iRow] = 1.0;
			}
			// perform gauss jordan elimination to see if the
			// identiry matrix appears on the LHS, and the RHS
			// is the inverse.
			//
			// 1 0 0 | J K L
			// 0 1 0 | M N O
			// 0 0 1 | P Q R
			//
			matResult = matResult.GaussJordonElimination();
			bool bRowHasAllZeros = true;
			// check for rows of zeros in the RHS
			for (int iRow = 0; iRow < iRows; iRow++)
			{
				bRowHasAllZeros = true;  // assume true until proven otherwise
				for (int iCol = 0; iCol < iColumns; iCol++)
				{
					if (matResult.arrMatrix[iRow, iColumns+iRow] != 0)
					{
						bRowHasAllZeros = false;
						break;
					}
				}
				// if there is a row with all zeros, then the matrix is not
				// invertible
				if (bRowHasAllZeros)
					break;
			}
			if (!bRowHasAllZeros)
			{
				Matrix matTemp= matResult;
				matResult = new Matrix(iRows, iColumns);
				for (int iRow = 0; iRow < iRows; iRow++)
				{
					for (int iCol = 0; iCol < iColumns; iCol++)
					{
						matResult.arrMatrix[iRow, iCol] = 
							matTemp.arrMatrix[iRow, iColumns+iCol];
					}
				}
			}
			else
			{
				matResult.Resize(0,0);
			}
			return matResult;
		}

		/// <summary>
		/// Makes the matrix an Identity matrix IF the matrix is
		/// in fact square.
		/// </summary>
		public void LoadIdentity()
		{
			if (arrMatrix == null)
				return;
			if (iRows != iColumns)
				throw new Exception("Cannot set as identity matrix - matrix is not square.");
			for (int iRow = 0; iRow < iRows; iRow++)
			{
				for (int iCol = 0; iCol < iColumns; iCol++)
				{
					if (iRow == iCol)
						arrMatrix[iRow, iCol] = 1;
					else
						arrMatrix[iRow, iCol] = 0;
				}
			}
		}

		/// <summary>
		/// Multiplies two matricies together and returns the result, 
		/// provided that the second matrix has the same number of rows 
		/// as the first matrix has columns. An exception is thrown 
		//  if the matricies are of the wrong size.
		/// </summary>
		public Matrix Multiply(Matrix mat)
		{
			if (arrMatrix == null)
				throw new NullReferenceException("Error in Matrix.Multiply." +
					"Provided matrix is null.");
			if (arrMatrix == null)
				throw new Exception("Error in Matrix.Multiply. " + 
					"This matrix is a null matrix.");
			if (mat.arrMatrix == null)
				throw new Exception("Error in Matrix.Multiply. " + 
					"Provided matrix is a null matrix.");
			if (iColumns != mat.iRows)
				throw new Exception("Error in Matrix.Multiply. " + 
					"The number of columns in this matrix is not equal " +
					"to the number of rows in the provided matrix.");
			// create a new result matrix
			int iResultRows = iRows;
			int iResultCols = mat.iColumns;
			Matrix matResult = new Matrix(iResultRows, iResultCols, true);
			for (int iResultRow = 0; iResultRow < iResultRows; iResultRow++)
			{
				for (int iResultCol = 0; iResultCol < iResultCols; iResultCol++)
				{
					for (int iCol = 0; iCol < iColumns; iCol++)
					{
						matResult.arrMatrix[iResultRow, iResultCol] += 
							arrMatrix[iResultRow, iCol] *
							mat.arrMatrix[iCol, iResultCol];
					}
				}
			}
			// return the results
			return matResult;
		}

		/// <summary>
		/// Multiplies a vector by the matrix in order to produce a new vector.
		/// An exception is thrown if the vector has a size that is not equal
		/// to then number of rows in the matrix, or if the matrix is null or
		/// vector is null. Equivalent to VECT * MAT = VECT
		/// </summary>
		public Vector Multiply(Vector vect)
		{
			if (vect == null)
				throw new NullReferenceException("Error in Matrix.Multiply." +
					"Provided vector object is null.");
			if (arrMatrix == null)
				throw new Exception("Error in Matrix.Multiply. " + 
					"This matrix is a null matrix.");
			if (vect.IsNull)
				throw new Exception("Error in Matrix.Multiply. " + 
					"Provided vector is a null vector.");
			int iRows = Rows;
			int iColumns = Columns;
			if (iRows != vect.Size)
				throw new Exception("Error in Matrix.Multiply. " + 
					"The number vector is not the same size as the number " +
					"of rows in the matrix.");
			// convert the vector to a matrix
			Matrix matResult = new Matrix(vect);
			matResult = matResult.Multiply(this);
			// create a new matrix that is of the same # or rows as the result
			Vector vectResult = new Vector(iColumns);
			for (int iIndex = 0; iIndex < iColumns; iIndex++)
			{
				vectResult[iIndex] = matResult[0, iIndex];
			}
			return vectResult;
		}

		/// <summary>
		/// Sets the dimensions of the matrix to those specified. Note
		/// that if either the row or column are <= 0, then the matrix
		/// becomes null. All data that is in the current matrix will
		/// be copied to the new matrix if possible.
		/// </summary>
		public void Resize(int iRows, int iColumns)
		{
			Resize(iRows, iColumns, false);
		}

		/// <summary>
		/// Sets the dimensions of the matrix to those specified. Note
		/// that if either the row or column are <= 0, then the matrix
		/// becomes null. All data that is in the current matrix will
		/// be copied to the new matrix if possible, unless bClearData
		/// is set to true.
		/// </summary>
		public void Resize(int iRows, int iColumns, bool bClearData)
		{
			if (iRows < 1 || iColumns < 1)
			{
				arrMatrix = null;
				return;
			}
			int iOldMatrixRows = this.iRows;
			int iOldMatrixColumns = this.iColumns;
			this.iRows = iRows;
			this.iColumns = iColumns;
			// allocate an array of rows
			double[,] arrNewMatrix = new double[iRows, iColumns];
			// in each row, allocate the array of columns within the row
			for (int iRow = 0; iRow < iRows; iRow++)
			{
				// transfer data from
				// the old matrix if possible and requested
				for (int iCol = 0; iCol < iColumns; iCol++)
				{
					if (arrMatrix != null && bClearData == false &&
						iRow < iOldMatrixRows && iCol < iOldMatrixColumns)
					{
						arrNewMatrix[iRow, iCol] = arrMatrix[iRow, iCol];
					}
					else
					{
						arrNewMatrix[iRow, iCol] = 0;
					}
				}
			}
			// use the new, up to date matrix
			arrMatrix = arrNewMatrix;
		}

		/// <summary>
		///Applies a scalar value to the matrix, which is the same as
		/// multiplying every value of the matrix by the scalar value
		/// provided.
		/// </summary>
		public void Scalar(double dScalar)
		{
			if (arrMatrix == null)
				return;
			int iRows = Rows;
			int iColumns = Columns;
			for (int iRow = 0; iRow < iRows; iRow++)
			{
				for (int iCol = 0; iCol < iColumns; iCol++)
				{
					arrMatrix[iRow, iCol] *= dScalar;
				}
			}
		}

		/// <summary>
		/// Outputs the context of teh matrix in an n x m dimensional
		/// spreadsheet.
		/// </summary>
		public override String ToString()
		{
			if (arrMatrix == null)
				return null;
			StringBuilder sbOutput = new StringBuilder();
			sbOutput.AppendLine("Matrix " + iRows + "x" + iColumns);
			for (int iRow = 0; iRow < iRows; iRow++)
			{
				for (int iCol = 0; iCol < iColumns; iCol++)
				{
					sbOutput.Append(String.Format("{0:G4}\t", arrMatrix[iRow, iCol]));
				}
				sbOutput.AppendLine();
			}
			return sbOutput.ToString();
		}

		/// <summary>
		/// Converts the matrix to the transpose of itself (i.e. reorients 
		/// and switches the contents of the row and column vectors).
		/// </summary>
		public Matrix Transpose()
		{
			Matrix matTransposed = new Matrix(iColumns, iRows, true);
			// copy the contents of the old into the new (iterating over the old matrix)
			for (int iRow = 0; iRow < iRows; iRow++)
			{
				for (int iCol = 0; iCol < iColumns; iCol++)
				{
					matTransposed.arrMatrix[iCol, iRow] = arrMatrix[iRow, iCol];
				}
			}
			// set the current matrix to the new one
			return matTransposed;
		}

		/// <summary>
		/// Zeros the matrix by setting all entries to zero. Same as Clear().
		/// </summary>
		public void Zero()
		{
			Clear();
		}

		/// <summary>
		/// Creates a new 2D array of doubles, without the overhead of initialization
		/// in Resize().
		/// </summary>
		protected double[,] Create2DArray(int iRows, int iColumns)
		{
			if (iRows < 0 || iColumns < 0)
				return null;
			double[,] arrNewMatrix = new double[iRows, iColumns];
			this.iRows = iRows;
			this.iColumns = iColumns;
			return arrNewMatrix;
		}
	};	// class Matrix

	public class MatrixMath
	{
		public static Vector CrossProduct(Vector vectA, Vector vectB)
		{
			if (vectA == null || vectB == null)
				throw new NullReferenceException();
			if (vectA.Size == 0 || vectB.Size == 0)
				throw new Exception("Error in MatrixMath.CrossProduct. Neither vector can be a null vector.");
			if (vectA.Size != 3 || vectB.Size != 3)
				throw new Exception("Error in MatrixMath.CrossProduct. Both vectors must be of size 3.");
			return new Vector(
				vectA[1] * vectB[2] - vectA[2] * vectB[1],
				vectA[2] * vectB[0] - vectA[0] * vectB[2],
				vectA[0] * vectB[1] - vectA[1] * vectB[0]);
		}

	}; // class MatrixMath

	public class UnitTesting
	{
		private Matrix matA;	   // 2x3 matrix
		private Matrix matB;      // 3x4 matrix
		private Matrix matC;      // 3x3 matrix
		private Matrix matNull;   // zero size matrix
		//private Matrix matZero;   // 3x3 matrix of zeros
		private Vector vectA;	   // vector of size 3
		private Vector vectB;	   // vector of size 3

		public UnitTesting()
		{
			matA = null;
			matB = null;
			matC = null;
			matNull = null;
			//matZero = null;
		}

		public void AllTests()
		{
			AllMatrixTests();
			AllVectorTests();
		}

		public void AllMatrixTests()
		{
			TestMatrixMultiply();
			TestMatrixTranspose();
			TestMatrixGaussJordanElimination();
			TestMatrixInverse();
		}

		public void AllVectorTests()
		{
			TestDotProduct();
			TestDotProductAngle();
			TestVectorMultiplication();
		}

		/// <summary>
		/// Creates matrix A:
		///
		///  1 2 4
		///  2 6 0
		///
		/// </summary>
		public void SetupMatA()
		{
			matA = new Matrix(2,3);
			matA[0, 0] = 1;
			matA[0, 1] = 2;
			matA[0, 2] = 4;

			matA[1, 0] = 2;
			matA[1, 1] = 6;
			matA[1, 2] = 0;
		}

		/// <summary>
		/// Creates matrix B:
		///
		///  4  1  4  3
		///  0 -1  3  1
		///  2  7  5  2
		/// </summary>
		public void SetupMatB()
		{
			matB = new Matrix(3,4);
			matB[0, 0] = 4;
			matB[0, 1] = 1;
			matB[0, 2] = 4;
			matB[0, 3] = 3;
			matB[1, 0] = 0;
			matB[1, 1] = -1;
			matB[1, 2] = 3;
			matB[1, 3] = 1;
			matB[2, 0] = 2;
			matB[2, 1] = 7;
			matB[2, 2] = 5;
			matB[2, 3] = 2;
		}

		/// <summary>
		/// Creates matrix C:
		///
		///  1  2  3
		///  2  5  3
		///  1  0  8
		/// </summary>
		public void SetupMatC()
		{
			matC = new Matrix(3,3);
			matC[0, 0] = 1;
			matC[0, 1] = 2;
			matC[0, 2] = 3;
			matC[1, 0] = 2;
			matC[1, 1] = 5;
			matC[1, 2] = 3;
			matC[2, 0] = 1;
			matC[2, 1] = 0;
			matC[2, 2] = 8;
		}


		/// <summary>
		/// Creates matrix Null, which a matrix of size 0.
		/// </summary>
		public void SetupMatNull()
		{
			matNull = new Matrix(0,0);
		}

		/// <summary>
		/// Creates a 3x3 zero matrix.
		/// </summary>
		public void SetupMatZero()
		{
			matC = new Matrix(3,3);
			matC[0, 0] = 0;
			matC[0, 1] = 0;
			matC[0, 2] = 0;
			matC[1, 0] = 0;
			matC[1, 1] = 0;
			matC[1, 2] = 0;
			matC[2, 0] = 0;
			matC[2, 1] = 0;
			matC[2, 2] = 0;
		}

		public void SetupVectA()
		{
			vectA = new Vector(3);
			vectA[0] = 2;
			vectA[1] = -1;
			vectA[2] = 1;
		}

		public void SetupVectB()
		{
			vectB = new Vector(3);
			vectB[0] = 1;
			vectB[1] = 1;
			vectB[2] = 2;
		}

		/// <summary>
		/// Tests the guass jordan elimination function of a matrix.
		/// </summary>
		public void TestMatrixGaussJordanElimination()
		{
			Debug.WriteLine("\n************  Test Gauss Jordon Elimination ****************");
			SetupMatB();
			Debug.WriteLine(matB.ToString());
			matC = matB.GaussJordonElimination();
			Debug.WriteLine("\nReduced Row Echelon Format is:");
			Debug.WriteLine(matC.ToString());
			Debug.WriteLine("\nReduced Row Echelon Format SHOULD be:");
			Debug.WriteLine("1       0       0       0.4556");
			Debug.WriteLine("0       1       0       -0.06667");
			Debug.WriteLine("0       0       1       0.3111");
		}

		/// <summary>
		/// Test matrix inversion.
		/// </summary>
		public void TestMatrixInverse()
		{
			Debug.WriteLine("\n************  Test Matrix Inverse ****************");
			SetupMatC();
			Debug.WriteLine(matC.ToString());
			Debug.WriteLine("\nC^-1 = \n");
			matC = matC.Inverse();
			Debug.WriteLine(matC.ToString());
			Debug.WriteLine("\nC^-1 SHOULD be:");
			Debug.WriteLine("-40    16       9");
			Debug.WriteLine("13     -5      -3");
			Debug.WriteLine("5      -2      -1");
		}

		/// <summary>
		/// Tests to see if two matricies can be properly mutlipled together.
		/// </summary>
		void TestMatrixMultiply()
		{
			Debug.WriteLine("\n************  Test Matrix Multiply  ****************");
			SetupMatA();
			Debug.WriteLine(matA.ToString());
			SetupMatB();
			Debug.WriteLine(matB.ToString());

			Debug.WriteLine("\nAB = C\n");
			matC = matA.Multiply(matB);
			Debug.WriteLine(matC.ToString());
			Debug.WriteLine("\nC SHOULD be:");
			Debug.WriteLine("12      27      30      13");
			
			Debug.WriteLine("8       -4      26      12");
		}

		/// <summary>
		/// Tests to see if the matrix can be properly transposed.
		/// </summary>
		void TestMatrixTranspose()
		{
			Debug.WriteLine("\n************  Test Matrix Transpose ****************");
			SetupMatB();
			Debug.WriteLine(matB.ToString());
			matC = matB.Transpose();
			Debug.WriteLine("\nTranspose is: \n");
			Debug.WriteLine(matC.ToString());
			Debug.WriteLine("\nC SHOULD be:");
			Debug.WriteLine("4       0       2");
			Debug.WriteLine("1       -1      7");
			Debug.WriteLine("4       3       5");
			Debug.WriteLine("3       1       2");
		}

		void TestDotProduct()
		{
			Debug.WriteLine("\n************  Test Vector Dot Product ****************");
			SetupVectA();
			Debug.WriteLine(vectA.ToString());
			SetupVectB();
			Debug.WriteLine(vectB.ToString());
			double d = vectA.DotProduct(vectB);
			Debug.WriteLine("Dot product is: " + d + "\n");
			Debug.WriteLine("Dot product SHOULD BE: 3\n");
		}

		void TestDotProductAngle()
		{
			Debug.WriteLine("\n************  Test Dot Product Angle ****************");
			SetupVectA();
			Debug.WriteLine(vectA.ToString());
			SetupVectB();
			Debug.WriteLine(vectB.ToString());
			double d = vectA.Angle(vectB);
			Debug.WriteLine("Dot product angle is: " + d + "\n");
			Debug.WriteLine("Dot product angle SHOULD BE: 60\n");
		}

		/// <summary>
		/// Tests to see if a vector can be multipled by a matrix
		/// </summary>
		void TestVectorMultiplication()
		{

		}


	};
}
