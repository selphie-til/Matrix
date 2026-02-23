//  Matrix.cpp

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <random>
#include <complex>

#include "matrix.hpp"

/**
 * @brief デストラクタ
 *
 * Matrix オブジェクトがスコープを抜ける際に呼ばれるデストラクタです。
 * Matrix オブジェクト内部で管理されている top_ メンバ変数への参照を null に設定します。
 *
 * @tparam T 行列内の各要素の型を定義するテンプレートパラメータ
 */
template<typename T>
Matrix<T>::~Matrix() {
    this->top_ = nullptr;
}

/**
 * @brief (ti,tj)がタイルの最後の場合には、行の端数を返し、それ以外の場合にはタイルの行数を返す。
 */
template<typename T>
uint32_t Matrix<T>::mb(const uint32_t &ti, const uint32_t &tj) const {
    assert( ti < this->p_);
    assert( tj < this->q_);
    return ( ( (this->m_) % (this->mb_) == 0) || (ti != ( (this->p_) - 1) ) ) ? (this->mb_) : (this->m_) % (this->mb_);
}

/**
 * @brief (ti,tj)がタイルの最後の場合には、列の端数を返し、それ以外の場合にはタイルの列数を返す。
 */
template<typename T>
uint32_t Matrix<T>::nb(const uint32_t &ti, const uint32_t &tj) const {
    assert( ti < this->p_);
    assert( tj < this->q_);
    return ( ( (this->n_) % (this->nb_) == 0) || (tj != ( (this->q_) - 1) ) ) ? (this->nb_) : (this->n_) % (this->nb_);
}

/**
 * @brief マトリクスの要素にランダムな数を割り当てます。
 *
 * この関数は標準の乱数生成器を使用してマトリクスのすべての要素にランダムな実数を割り当てます。
 *
 * @tparam T マトリクスの要素の型
 *
 * 乱数の範囲は0.0から1.0までです。
 */
template<typename T>
void Matrix<T>::gen_rnd_elm() {
    // 乱数生成器
    std::random_device rdev;
    // ランダムなシードの設定
    std::mt19937 engine(rdev());
    std::uniform_real_distribution<T> dist(0.0, 1.0);
    for (uint64_t i = 0; i < m_*n_; ++i) {
        ((this->top_).get())[i] = dist(engine);
    }
}

/**
 * @brief 行列の特定の要素を指すポインタを返します。
 *
 * @param ti タイルの行インデックス。
 * @param tj タイルの列インデックス。
 * @return 行列の(ti, tj)位置にある要素へのポインタ。
 *
 * tiとtjは行列の行と列を示しており、0から始まるインデックス番号として解釈されます。関数は、指定されたインデックス位置の行列要素へのポインタを返します。
 *
 * tiとtjが行列の有効な範囲内ではない場合、アサーションでエラーが引き起こされます。
 */
template<typename T>
T* Matrix<T>::elm(const uint32_t &ti, const uint32_t &tj) const {
    assert(ti < (this->p_));
    assert(tj < (this->q_));

    uint64_t pos = Matrix<T>::convertTileToArray(ti, tj, 0, 0);

    return &(top_.get()[pos]);
}

/**
 * @brief Matrixクラス内の特定の要素へのポインタを返します。
 *
 * @tparam T マトリックスのデータ型
 * @param ti マトリックスの行方向のタイルインデックス
 * @param tj マトリックスの列方向のタイルインデックス
 * @param i 行インデックス
 * @param j 列インデックス
 * これらのパラメータは、マトリックス内の特定の要素を指定します。
 *
 * @return T* 特定の要素へのポインタを返します。
 *
 * @details
 * 最初に、関数は行と列のインデックスが適切な範囲内にあることを検証します。
 * 適切な位置が計算され、最終的にデータの要素へのポインタが返されます。
 * タイル内のインデックスを指定することにより、行列が大きい場合でも効率的にデータにアクセスできます。
 */
template<typename T>
T* Matrix<T>::elm(const uint32_t &ti, const uint32_t &tj,
                  const uint32_t &i, const uint32_t &j) const {
    assert(i < (ti == (this->p_ - 1) ? ((this->m_) % (this->mb_) == 0) ? this->mb_ : (this->m_) % (this->mb_) : this->mb_));
    assert(j < (tj == (this->q_ - 1) ? ((this->n_) % (this->nb_) == 0) ? this->nb_ : (this->n_) % (this->nb_) : this->nb_));

    assert(ti < this->p_);
    assert(tj < this->q_);

    uint64_t pos = Matrix<T>::convertTileToArray( ti, tj, i, j);

    return &(top_.get()[pos]);
}

/**
 * @brief マトリックスデータを指定のファイルに書き出すメソッド
 *
 * @tparam T マトリックス要素の型
 * @param fname [in] 書き出し先のファイル名
 *
 * @details
 * fnameに指定されたファイルにマトリックスデータを書き出します。各要素の値は空白で区切ります。
 * ファイルオープンに失敗した場合は、エラーメッセージを出力して終了します。
 *
 * 内部では精度5で出力を行い、行と列の順番に従って、マトリックスの全ての要素を順番に書き出します。
 * ファイルへの書き出しが終わったら、ファイルをクローズします。
 */
template<typename T>
void Matrix<T>::file_out(const char *fname) {

    std::ofstream matf(fname);
    if (!matf) {
        std::cerr << "Unable to open " << fname << std::endl;
        exit(1);
    }

    matf.precision(5);

    for (uint32_t i = 0; i < this->m_; i++) {
        for (uint32_t j = 0; j < this->n_; j++) {
            const uint32_t ti = i / this->mb_;
            const uint32_t tj = j / this->nb_;
            const uint32_t tp = i % this->mb_;
            const uint32_t tq = j % this->nb_;

            uint64_t pos = Matrix<T>::convertTileToArray( ti, tj, tp, tq);

            matf << (this->top_)[pos] << " ";
        }
        matf << std::endl;
    }
    matf.close();
}

/**
 * @brief 代入演算子のオーバーロード
 *
 * この演算子オーバーロードは、引数として渡された 'Matrix<T>' 型のオブジェクトの値を、このオブジェクトにコピーします。
 * 以下のメンバ変数がコピーされます: m_, n_, mb_, nb_, p_, q_, top_.
 *
 * @param M コピー元の 'Matrix<T>' 型のオブジェクト
 * @return コピー後の自身のオブジェクトの参照
 * @throw assert アサーションに失敗した場合にスローされます。(this->m_ != M.m_ または this->n_ != M.n_ の場合)
 */
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &M){
    assert(this->m_ == M.m_);
    assert(this->n_ == M.n_);

    this->m_ = M.m_;
    this->n_ = M.n_;
    this->mb_ = M.mb_;
    this->nb_ = M.nb_;
    this->p_ = M.p_;
    this->q_ = M.q_;

    for (uint64_t i = 0; i < (this->m_)*(this->n_); i++)
        (this->top_).get()[i] = (M.top_).get()[i];

    return *this;
}

/**
 * @brief      Matrixクラスの要素同士を加算するメソッド。
 *
 * @tparam     T     Matrixが持つデータ型。
 * @param[in]  M     加算対象となるMatrixオブジェクトの参照。
 *
 * @return     二つのMatrixの各要素を加算した新たなMatrixオブジェクトを返します。
 *
 * @details    このメソッドは自身のMatrixと引数で渡されたMatrixの各要素を加算します。ただし、この二つのMatrixの次元は同じでなければいけません。
 *             それぞれの要素を加算した結果を新たなMatrixオブジェクトに格納し、そのオブジェクトを返します。
 */
template<typename T>
Matrix<T> Matrix<T>::addElements(const Matrix<T> &M) const {
    assert(this->m_ == M.m_);
    assert(this->n_ == M.n_);

    Matrix N(this->m_, this->n_, this->mb_, this->nb_);
    uint64_t totalElements = (this->m_) * (this->n_);
    for( uint64_t i=0; i<totalElements; ++i){
        N[i] = (*this)[i] + M[i];
    }

    return N;
}

/**
 * @brief マトリックスの加算を行う演算子オーバーロード
 *
 * このオーバーロードされた+演算子は、指定されたMatrixオブジェクト（M）と現在のMatrixオブジェクト間で加算を行います。
 * ただし、このメソッドを呼び出す前に、両方のMatrixが同じ次元を持つことが確認されていることが期待されています（そうでない場合はアサートに失敗します）。
 *
 * @param M 現在のMatrixオブジェクトに加算するMatrixオブジェクト
 * @return 2つのMatrixオブジェクトの合計を持つ新しいMatrixオブジェクト
 * @throw assertion failed 両方のMatrixオブジェクトが同じ次元でない場合、アサートに失敗します。
 *
 * @tparam T Matrix要素の型
 */
template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &M) const {
    assert(m_ == M.m_);
    assert(n_ == M.n_);

    return this->addElements(M);
}

/**
 * @brief この関数は、引数として与えられた行列の要素を対応する要素に減算する新しい行列を作成します。
 *
 * @tparam T 行列の要素の型
 * @param M 減算を行う対象の行列。今の行列と同じ次元である必要があります。
 * @return Matrix<T> 減算結果をもつ新しい行列を返します。
 *
 * @pre この関数は、対象の行列が同じ次元であることを必須とします。
 * @post 減算結果をもつ新しい行列が返されます。
 *
 * @note もし渡された行列の次元が異なる場合、アサートエラーを返します。
 */
template<typename T>
Matrix<T> Matrix<T>::minusElements(const Matrix<T> &M) const{
    assert(this->m_ == M.m_);
    assert(this->n_ == M.n_);

    Matrix N( this->m_, this->n_, this->mb_, this->nb_);
    uint64_t totalElements = (this->m_) * (this->n_);
    for( uint64_t i=0; i<totalElements; ++i){
        N[i] = (*this)[i] - M[i];
    }

    return N;
}

/**
 * @brief 行列の減算を行います。
 *
 * このメソッドは、現在のMatrix<T>オブジェクトから別のMatrix<T>オブジェクトを減算します。
 * 減算は要素ごとに行われ、新しいMatrix<T>オブジェクトが返されます。
 * 減算を行う前に、行列の大きさが同じであることが確認されます。
 *
 * @tparam T 行列の要素の型
 * @param M 減算するMatrix<T>オブジェクト
 * @return 新しいMatrix<T>オブジェクト、計算結果
 */
template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &M) const {
    assert(m_ == M.m_);
    assert(n_ == M.n_);

    return this->minusElements(M);
}

/**
 * @brief 与えられた行列とこの行列が等しいかどうかをチェックするためのオーバーロードされた等価演算子
 *
 * @tparam T 行列に格納される要素の型
 * @param M 比較対象のMatrixオブジェクトへの参照
 * @return bool この行列と与えられた行列が等しい場合はtrue、そうでない場合はfalseを返します。
 *
 * この演算子はまず、両方の行列が同じ形状であることを確認します（つまり、行数と列数が一致します）。
 * 形状が一致しない場合、演算子はfalseを即座に返します。
 * その後、行列のすべての要素を順番に比較し、一致しない要素がある場合、演算子はfalseを返し、
 * すべての要素が一致する場合、演算子はtrueを返します。
 */
template<typename T>
bool Matrix<T>::operator==(const Matrix<T> &M) {
    if(m_ != M.m_ || n_ != M.n_) {
        return false;
    }
    for (uint64_t i = 0; i < m_*n_; i++) {
        if (this->top_[i] != M.top_[i]){
            return false;
        }
    }
    return true;
}

/**
 * @brief オーバーロードした[]演算子
 *
 * @tparam T マトリックスの要素の型
 * @param i マトリックス要素の一次元配列表現におけるインデックス
 * @return リファレンス型T : マトリックスのi番目の要素へのリファレンス
 *
 * @details
 * この関数は、マトリックスの要素へのアクセスを提供するために[]演算子をオーバーロードします。
 * 一次元配列表現におけるインデックスを引数に取り、該当するマトリックス要素へのリファレンスを返します。
 * 入力パラメータ i は 0 以上 (m_)*(n_) 未満であることをアサートで確認します。
 */
template<typename T>
T &Matrix<T>::operator[](const uint64_t &i) const {
    assert(i < (this->m_)*(this->n_));

    return top_[i];
}
/**
 * @brief オーバーロードされた()演算子を使ってマトリクス内の特定の要素を取得します。
 *
 * @tparam T マトリクスが保持するデータの型
 *
 * @param i マトリクス内の行の位置
 * @param j マトリクス内の列の位置
 * @return T& 指定された位置(i, j)のマトリクスの要素への参照
 *
 * 当該要素の位置がマトリクスの範囲外であると、assertでプログラムが終了します。
 */
template<typename T>
T &Matrix<T>::operator()(const uint32_t &i, const uint32_t &j) const {
    assert(i < (this->m_));
    assert(j < (this->n_));

    const uint32_t ti = i / this->mb_;
    const uint32_t tj = j / this->nb_;
    const uint32_t tp = i % this->mb_;
    const uint32_t tq = j % this->nb_;

    uint64_t pos = Matrix<T>::convertTileToArray( ti, tj, tp, tq);

    return top_[ pos ];
}

/**
 * @brief すべてのマトリックス要素を0に設定します。
 *
 * この関数は内部の`top_`メンバで指定された配列に対してstd::fillを使用し、
 * `m_` * `n_` 個の要素すべてを0に設定します。
 *
 * @tparam T マトリックス要素のタイプ
 */
template<typename T>
void Matrix<T>::zero() {

    std::fill(top_.get(), top_.get() + (this->m_) + (this->n_), 0);

}

/**
 * @brief タイル行列から配列への変換を行う
 *
 * タイル行列の座標 (ti, tj) および、タイル内の座標 (i, j) から、配列のインデックスを計算します。
 * また、タイル行列の配置（Column-Major/Row-Major）および、タイル内の配置（Column-Major/Row-Major）に基づいた順序付けを制御します。
 *
 * @tparam T 値の型
 * @param ti タイル行列の行番号
 * @param tj タイル行列の列番号
 * @param i タイル内の行番号
 * @param j タイル内の列番号
 * @return 生成された配列のインデックス
 */
template<typename T>
uint64_t Matrix<T>::convertTileToArray(const uint32_t &ti, const uint32_t &tj, const uint32_t &i, const uint32_t &j) const{
    uint64_t index = {0};
    switch ( static_cast<std::underlying_type<Ordering>::type>(ordering_) ) {
        case static_cast<std::underlying_type<Ordering>::type>(Ordering::TileMatrixColumnMajorTileColumnMajor):
            index += static_cast<uint64_t>(ti) * static_cast<uint64_t>(this->mb_) * static_cast<uint64_t>(Matrix<T>::nb(ti, tj));
            index += static_cast<uint64_t>(tj) * static_cast<uint64_t>(this->m_) * static_cast<uint64_t>(this->nb_);
            index += static_cast<uint64_t>(j) * static_cast<uint64_t>(Matrix<T>::mb(ti, tj));
            index += static_cast<uint64_t>(i);
            break;
        case static_cast<std::underlying_type<Ordering>::type>(Ordering::TileMatrixRowMajorTileColumnMajor):
            index += static_cast<uint64_t>(ti) * static_cast<uint64_t>(this->mb_) * static_cast<uint64_t>(this->n_);
            index += static_cast<uint64_t>(tj) * static_cast<uint64_t>(Matrix<T>::mb(ti, tj)) * static_cast<uint64_t>(this->nb_);
            index += static_cast<uint64_t>(j) * static_cast<uint64_t>(Matrix<T>::mb(ti, tj));
            index += static_cast<uint64_t>(i);
            break;
        case static_cast<std::underlying_type<Ordering>::type>(Ordering::TileMatrixColumnMajorTileRowMajor):
            index += static_cast<uint64_t>(ti) * static_cast<uint64_t>(this->mb_) * static_cast<uint64_t>(Matrix<T>::nb(ti, tj));
            index += static_cast<uint64_t>(tj) * static_cast<uint64_t>(this->m_) * static_cast<uint64_t>(this->nb_);
            index += static_cast<uint64_t>(i) * static_cast<uint64_t>(Matrix<T>::nb(ti, tj));
            index += static_cast<uint64_t>(j);
            break;
        case static_cast<std::underlying_type<Ordering>::type>(Ordering::TileMatrixRowMajorTileRowMajor):
            index += static_cast<uint64_t>(ti) * static_cast<uint64_t>(this->mb_) * static_cast<uint64_t>(this->n_);
            index += static_cast<uint64_t>(tj) * static_cast<uint64_t>(Matrix<T>::mb(ti, tj)) * static_cast<uint64_t>(this->nb_);
            index += static_cast<uint64_t>(i) * static_cast<uint64_t>(Matrix<T>::nb(ti, tj));
            index += static_cast<uint64_t>(j);
            break;
        default:
            break;
    }

    assert( index < static_cast<uint64_t>((this->m_)) * static_cast<uint64_t>((this->n_)));

    return index;
}

/**
 * @brief set_ts メソッドは、新しいタイルサイズを設定します。
 *
 * @details この関数は、指定したタイルサイズで新しい行列を作成し、元の行列のデータを新しい行列にコピーします。
 * 新しい行列が作成された後、元の行列は新しい行列に置き換えられます。
 * タイルサイズが0未満の場合はアサートで失敗します。
 *
 * @param ts 新しいタイルサイズ。
 */
template<typename T>
void Matrix<T>::set_ts(uint32_t ts) {
    assert( ts < m_ );
    assert( ts < n_ );

    Matrix new_matrix( m_, n_, ts, ordering_);

    for( uint32_t i=0; i<m_; ++i){
        for(uint32_t j=0; j<n_; ++j){
            new_matrix(i,j) = this->operator()(i,j);
        }
    }

    *this = std::move(new_matrix);
}

/**
   * @brief Matrix<T>::set_ts() メソッドは、Matrix<T> オブジェクトのタイルサイズを設定します。
   *
   * @param mb    新しい行のタイルサイズ (1より大きい必要があります)
   * @param nb    新しい列のタイルサイズ (1より大きい必要があります)
   */
template<typename T>
void Matrix<T>::set_ts(uint32_t mb, uint32_t nb) {
    assert( mb <= m_);
    assert( nb <= n_);

    Matrix new_matrix( m_, n_, mb, nb, ordering_);

    for( uint32_t i=0; i<m_; ++i){
        for(uint32_t j=0; j<n_; ++j){
            new_matrix(i,j) = this->operator()(i,j);
        }
    }

    *this = std::move(new_matrix);
}


//! フロート値用のMatrix特化クラス
/*! このクラスはMatrixクラスのfloat特化版です。具体的な演算の実装を行っています。 */
template class Matrix<float>;

//! ダブル値用のMatrix特化クラス
/*! このクラスはMatrixクラスのdouble特化版です。具体的な演算の実装を行っています。 */
template class Matrix<double>;
