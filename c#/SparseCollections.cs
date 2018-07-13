// CREDIT TO:
// http://www.blackbeltcoder.com/Articles/algorithms/creating-a-sparse-matrix-in-net
// Jonathan Wood
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SparseCollections {
    public class SparseMatrixJonathanWood<T>
    {
        // Master dictionary hold rows of column dictionary
        protected Dictionary<int, Dictionary<int, T>> _rows;

        /// <summary>
        /// Constructs a SparseMatrix instance.
        /// </summary>
        public SparseMatrixJonathanWood()
        {
            _rows = new Dictionary<int, Dictionary<int, T>>();
        }

        /// <summary>
        /// Gets or sets the value at the specified matrix position.
        /// </summary>
        /// <param name="row">Matrix row</param>
        /// <param name="col">Matrix column</param>
        public T this[int row, int col]
        {
            get
            {
                return GetAt(row, col);
            }
            set
            {
                SetAt(row, col, value);
            }
        }

        /// <summary>
        /// Gets the value at the specified matrix position.
        /// </summary>
        /// <param name="row">Matrix row</param>
        /// <param name="col">Matrix column</param>
        /// <returns>Value at the specified position</returns>
        public T GetAt(int row, int col)
        {
            Dictionary<int, T> cols;
            if (_rows.TryGetValue(row, out cols))
            {
                T value = default(T);
                if (cols.TryGetValue(col, out value))
                    return value;
            }
            return default(T);
        }

        /// <summary>
        /// Sets the value at the specified matrix position.
        /// </summary>
        /// <param name="row">Matrix row</param>
        /// <param name="col">Matrix column</param>
        /// <param name="value">New value</param>
        public void SetAt(int row, int col, T value)
        {
            if (EqualityComparer<T>.Default.Equals(value, default(T)))
            {
                // Remove any existing object if value is default(T)
                RemoveAt(row, col);
            }
            else
            {
                // Set value
                Dictionary<int, T> cols;
                if (!_rows.TryGetValue(row, out cols))
                {
                    cols = new Dictionary<int, T>();
                    _rows.Add(row, cols);
                }
                cols[col] = value; }
        }

        /// <summary>
        /// Removes the value at the specified matrix position.
        /// </summary>
        /// <param name="row">Matrix row</param>
        /// <param name="col">Matrix column</param>
        public void RemoveAt(int row, int col)
        {
            Dictionary<int, T> cols;
            if (_rows.TryGetValue(row, out cols))
            {
                // Remove column from this row
                cols.Remove(col);
                // Remove entire row if empty
                if (cols.Count == 0)
                    _rows.Remove(row);
            }
        }

        /// <summary>
        /// Returns all items in the specified row.
        /// </summary>
        /// <param name="row">Matrix row</param>
        public IEnumerable<T> GetRowData(int row)
        {
            Dictionary<int, T> cols;
            if (_rows.TryGetValue(row, out cols))
            {
                foreach (KeyValuePair<int, T> pair in cols)
                {
                    yield return pair.Value;
                }
            }
        }

        /// <summary>
        /// Returns the number of items in the specified row.
        /// </summary>
        /// <param name="row">Matrix row</param>
        public int GetRowDataCount(int row)
        {
            Dictionary<int, T> cols;
            if (_rows.TryGetValue(row, out cols))
            {
                return cols.Count;
            }
            return 0;
        }

        /// <summary>
        /// Returns all items in the specified column.
        /// This method is less efficent than GetRowData().
        /// </summary>
        /// <param name="col">Matrix column</param>
        /// <returns></returns>
        public IEnumerable<T> GetColumnData(int col)
        {
            foreach (KeyValuePair<int, Dictionary<int, T>> rowdata in _rows)
            {
                T result;
                if (rowdata.Value.TryGetValue(col, out result))
                    yield return result;
            }
        }

        /// <summary>
        /// Returns the number of items in the specified column.
        /// This method is less efficent than GetRowDataCount().
        /// </summary>
        /// <param name="col">Matrix column</param>
        public int GetColumnDataCount(int col)
        {
            int result = 0;

            foreach (KeyValuePair<int, Dictionary<int, T>> cols in _rows)
            {
                if (cols.Value.ContainsKey(col))
                    result++;
            }
            return result;
        }
    }

    public class SparseMatrix<T> : SparseMatrixJonathanWood<T>, IEnumerable<T> {

        public SparseMatrix() : base() {}

        public T this[Tuple<int, int> ij] {
            get {
                return GetAt(ij.Item1, ij.Item2);
            }
            set {
                SetAt(ij.Item1, ij.Item2, value);
            }
        }

        public int Count() {
            int L=0;
            foreach (var row in _rows.Values) L += row.Count;
            return L;
        }

        /// <summary>
        /// Returns all row indices.
        /// </summary>
        public IEnumerable<int> GetRowIndices() {
            foreach (var row in _rows) {
                yield return row.Key;
            }
        }

        /// <summary>
        /// Returns all row indices.
        /// Slower than GetRowIndices()
        /// </summary>
        public IEnumerable<int> GetColumnIndices() {
            var columnIndices = new HashSet<int>();
            foreach (var ij in IndexPairs) columnIndices.Add(ij.Item2);
            foreach (int colIdx in columnIndices) {
                yield return colIdx;
            }
        }

        /// <summary>
        /// Returns all column indices in a specified row.
        /// </summary>
        /// <param name="rowIdx">Matrix row</param>
        public IEnumerable<int> GetColumnIndicesInRow(int rowIdx) {
            Dictionary<int, T> cols;
            if (_rows.TryGetValue(rowIdx, out cols)) {
                foreach (int colIdx in cols.Keys) {
                    yield return colIdx;
                }
            }
        }

        /// <summary>
        /// Returns all column indices in a specified row.
        /// Slower than GetColumnIndicesInRow(int)
        /// </summary>
        /// <param name="colIdx">Matrix column</param>
        public IEnumerable<int> GetRowIndicesInColumn(int colIdx) {
            foreach (KeyValuePair<int, Dictionary<int, T>> rowdata in _rows) {
                T result;
                if (rowdata.Value.TryGetValue(colIdx, out result))
                    yield return rowdata.Key;
            }
        }

        public IEnumerable<Tuple<int, int>> IndexPairs {
            get {
                foreach (var row in _rows) {
                    foreach (var col in row.Value) {
                        yield return new Tuple<int, int>(row.Key,col.Key);
                    }
                }
            }
        }



        public SparseMatrix<T> GetTranspose() {
            SparseMatrix<T> transpose = new SparseMatrix<T>();
            foreach (int rowIdx in _rows.Keys) {
                foreach (int colIdx in _rows[rowIdx].Keys) {
                    transpose[colIdx, rowIdx] = this[rowIdx, colIdx];
                }
            }
            return transpose;
        }
        public double Connectance() {
            var indices = new HashSet<int>();
            int numLinks = 0;
            foreach (var ij in IndexPairs) {
                indices.Add(ij.Item1);
                indices.Add(ij.Item2);
                numLinks++;
            }
            int numSpecies = indices.Count;
            return (double)numLinks / (numSpecies*numSpecies);
        }

        public string ToString(Func<T, string> Str) {
            var sb = new StringBuilder("{");
            foreach (int i in _rows.Keys.OrderBy(x=>x)) {
                sb.Append("{").Append(i.ToString()).Append(": {");
                foreach (int j in _rows[i].Keys.OrderBy(x=>x)) {
                    sb.Append(j.ToString()).Append(": ").Append(Str(GetAt(i,j))).Append(", ");
                }
                sb.Length -= 2;
                sb.Append("},\n");
            }
            sb.Length -= 2;
            sb.Append("}");
            return sb.ToString();
        }
        public int MaxIdx() {
            int max = 0;
            foreach (int i in _rows.Keys) {
                max = i > max ? i : max;
                foreach (int j in _rows[i].Keys) {
                    max = j > max ? j : max;
                }
            }
            return max;
        }
        public T[,] ToArray(int minIdx, int range) {
            if (range < 1) throw new System.InvalidOperationException("range < 1");
            var arr = new T[range,range];
            for (int i=0; i<range; i++) {
                for (int j=0; j<range; j++) {
                    arr[i+minIdx, j+minIdx] = GetAt(i,j);
                }
            }
            return arr;
        }

        public IEnumerator<T> GetEnumerator() {
            foreach (var row in _rows) {
                foreach (var col in row.Value) {
                    yield return col.Value;
                }
            }
        }
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() {
            return GetEnumerator();
        }
    }

    public class SparseVector<T> : IEnumerable<T> {
        Dictionary<int, T> _items;

        public SparseVector() {
            _items = new Dictionary<int, T>();
        }
        public SparseVector(IEnumerable<T> items) : this() {
            int i=0;
            foreach (T item in items) _items[i++] = item;
        }

        // <summary>
        // Gets or sets the value at the specified matrix position.
        // </summary>
        // <param name="row">Matrix row</param>
        // <param name="col">Matrix column</param>
        public T this[int idx] {
            get {
                return GetAt(idx);
            }
            set {
                SetAt(idx, value);
            }
        }

        // <summary>
        // Gets the value at the specified vector position.
        // </summary>
        // <param name="idx">Vector index</param>
        // <returns>Value at the specified position</returns>
        public T GetAt(int idx) {
            T value = default(T);
            if (_items.TryGetValue(idx, out value))
                return value;
            return default(T);
        }

        // <summary>
        // Sets the value at the specified vector position.
        // </summary>
        // <param name="idx">Vector Index</param>
        // <param name="value">New value</param>
        public void SetAt(int idx, T value) {
            if (EqualityComparer<T>.Default.Equals(value, default(T))) {
                // Remove any existing object if value is default(T)
                RemoveAt(idx);
            } else {
                // Set value
                _items[idx] = value;
            }
        }

        // <summary>
        // Removes the value at the specified matrix position.
        // </summary>
        // <param name="row">Matrix row</param>
        // <param name="col">Matrix column</param>
        public void RemoveAt(int idx) {
            _items.Remove(idx);
        }

        public IEnumerable<T> Data {
            get {
                foreach (T item in _items.Values) {
                    yield return item;
                }
            }
        }
        public IEnumerable<int> Indices {
            get {
                foreach (int i in _items.Keys) {
                    yield return i;
                }
            }
        }
        public int Count {
            get {
                return _items.Count;
            }
        }
        public IEnumerator<T> GetEnumerator() {
            return _items.Values.GetEnumerator();
        }
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() {
            return GetEnumerator();
        }
    }
}
