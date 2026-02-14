//  Matrix.hpp
#pragma once

#include <iostream>
#include <cstdint>
#include <memory>

/**
 * @enum Ordering
 * @brief タイル化された行列とその行列上のタイルのストレージ順序を表す列挙体。
 *
 * @var Ordering::TileMatrixRowMajorTileRowMajor
 * タイル行列自身が行方向優先であり、各タイル内の要素も行方向優先で格納されています。
 *
 * @var Ordering::TileMatrixRowMajorTileColumnMajor
 * タイル行列自身が行方向優先で、各タイル内の要素は列方向優先で格納されています。
 *
 * @var Ordering::TileMatrixColumnMajorTileRowMajor
 * タイル行列自身が列方向優先で、各タイル内の要素は行方向優先で格納されています。
 *
 * @var Ordering::TileMatrixColumnMajorTileColumnMajor
 * タイル行列自身も各タイル内も列方向優先で格納されています。
 */
enum class Ordering :unsigned short{
    TileMatrixRowMajorTileRowMajor,
    TileMatrixRowMajorTileColumnMajor,
    TileMatrixColumnMajorTileRowMajor,
    TileMatrixColumnMajorTileColumnMajor,
};

template<typename T>
class Matrix;

template<typename T>
std::ostream& operator<< ( std::ostream& os, Matrix<T> const &ma);

/**
 * @tparam T テンプレート変数。行列の要素の型。
 * @class Matrix
 * @brief 行列を表現するクラス。
 */
template <typename T>
class Matrix {

protected:
    /**
     * @brief 行列の行数。
     */
    uint32_t m_;  // number of lows of the matrix or lda
    /**
     * @brief 行列の列数。
     */
    uint32_t n_;  // number of columns of the matrix
    /**
     * @brief タイルの行数。
     */
    uint32_t mb_; // number of lows of the tile
    /**
     * @brief タイルの列数。
     */
    uint32_t nb_; // number of columns of the tile
    /**
     * @brief 行方向のタイル数。
     */
    uint32_t p_;  // number of low tiles
    /**
     * @brief 列方向のタイル数。
     */
    uint32_t q_;  // number of column tiles

    /**
     * @brief 行列の先頭への一意なポインタ。
     */
    std::unique_ptr<T[]> top_;     // pointer to the matrix

    /**
     * @brief オーダリングオブジェクト
     *
     * このメンバは、タイル行列の配置順序を管理します。
     * タイル行列、タイル内のrow-majorまたはcolumn-majorのいずれかになります。
     */
    Ordering ordering_;

public:
    /**
     * @brief デフォルトコンストラクタ。
     */
    Matrix():m_{0}, n_{0},
             mb_{0}, nb_{0},
             p_{0}, q_{0}, top_{nullptr},
             ordering_{Ordering::TileMatrixColumnMajorTileColumnMajor}
             {};
    /**
     * @brief 行列のサイズを指定するコンストラクタ。
     * @param m 行数
     * @param n 列数
     */
    Matrix(const uint32_t &m, const uint32_t &n,
           const Ordering ordering=Ordering::TileMatrixColumnMajorTileColumnMajor)
           :m_{m}, n_{n},
            mb_{m}, nb_{n},
            p_{1}, q_{1}, top_{std::make_unique<T[]>(static_cast<uint64_t>(m)*static_cast<uint64_t>(n))},
            ordering_{ordering} {};
    /**
     * @brief 行列のサイズとタイルのサイズを指定するコンストラクタ。
     * @param m 行数
     * @param n 列数
     * @param ts タイルのサイズ
     */
    Matrix(const uint32_t &m, const uint32_t &n, const uint32_t &ts,
           const Ordering ordering=Ordering::TileMatrixColumnMajorTileColumnMajor)
           :m_{m}, n_{n},
            mb_{ts}, nb_{ts},
            p_{(m % ts == 0) ? m / ts : m / ts + 1},
            q_{(n % ts == 0) ? n / ts : n / ts + 1},
            top_{std::make_unique<T[]>(static_cast<uint64_t>(m)*static_cast<uint64_t>(n))},
            ordering_{ordering} {};
    /**
     * @brief  行列のサイズとタイルのサイズを指定するコンストラクタ。
     * @param m 行数
     * @param n 列数
     * @param mb タイルの行数
     * @param nb タイルの列数
     */
    Matrix(const uint32_t &m, const uint32_t &n,
           const uint32_t &mb, const uint32_t &nb,
           const Ordering ordering=Ordering::TileMatrixColumnMajorTileColumnMajor):
            m_{m}, n_{n},
            mb_{mb}, nb_{nb},
            p_{(m % mb == 0) ? m / mb : m / mb + 1},
            q_{(n % nb == 0) ? n / nb : n / nb + 1},
            top_{std::make_unique<T[]>(static_cast<uint64_t>(m)*static_cast<uint64_t>(n))},
            ordering_{ordering} {};

    /**
     * @brief デストラクタ。
     */
    ~Matrix();

    /**
     * @brief コピーコンストラクタ。
     * @param M コピー元のMatrixオブジェクト
     */
    Matrix(const Matrix& M):m_{M.m_}, n_{M.n_}, mb_{M.mb_}, nb_{M.nb_},
                            p_{M.p_}, q_{M.q_}, top_{std::make_unique<T[]>(M.m_*M.n_)},
                            ordering_{M.ordering_}
    {
        if( M.top_ != nullptr) {
            std::copy((M.top_).get(), (M.top_).get() + M.m_ * M.n_, top_.get());
        }
    }

    /**
     * @brief ムーブコンストラクタ。
     * @param M ムーブ元のMatrixオブジェクト
     */
    Matrix(Matrix<T>&& M) noexcept: m_{M.m_}, n_{M.n_}, mb_{M.mb_}, nb_{M.nb_},
                                    p_{M.p_}, q_{M.q_}, top_{std::move(M.top_)},
                                    ordering_{M.ordering_}
    {
        M.m_ = 0;
        M.n_ = 0;
        M.mb_ = 0;
        M.nb_ = 0;
        M.p_ = 0;
        M.q_ = 0;
    }

    /**
     * @brief 行列の先頭へのポインタを取得するゲッター。
     */
    T* top() { return (this->top_).get(); }
    /**
     * @brief 行数を取得するゲッター。
     */
    [[nodiscard]] uint32_t m() const { return this->m_; }
    /**
     * @brief 列数を取得するゲッター。
     */
    [[nodiscard]] uint32_t n() const { return this->n_; }
    /**
     * @brief タイルの行数を取得するゲッター。
     */
    [[nodiscard]] uint32_t mb() const { return this->mb_; }
    /**
     * @brief (ti,tj)がタイルの最後の場合には、行の端数を返し、それ以外の場合にはタイルの行数を返す。
     */
    [[nodiscard]] uint32_t mb(const uint32_t &ti, const uint32_t &tj) const;
    /**
     * @brief タイルの列数を取得するゲッター。
     */
    [[nodiscard]] uint32_t nb() const { return this->nb_; }
    /**
     * @brief (ti,tj)がタイルの最後の場合には、列の端数を返し、それ以外の場合にはタイルの列数を返す。
     */
    [[nodiscard]] uint32_t nb(const uint32_t &ti, const uint32_t &tj) const;
    /**
     * @brief 行方向のタイル数を取得するゲッター。
     */
    [[nodiscard]] uint32_t p() const { return this->p_; }
    /**
     * @brief 列方向のタイル数を取得するゲッター。
     */
    [[nodiscard]] uint32_t q() const { return this->q_; }

    /**
     * @brief set_ts メソッドは、新しいタイルサイズを設定します。
     *
     * @details この関数は、指定したタイルサイズで新しい行列を作成し、元の行列のデータを新しい行列にコピーします。
     * 新しい行列が作成された後、元の行列は新しい行列に置き換えられます。
     * タイルサイズが0未満の場合はアサートで失敗します。
     *
     * @param ts 新しいタイルサイズ。
     */
    void set_ts(uint32_t ts);

    /**
     * @brief Matrix<T>::set_ts() メソッドは、Matrix<T> オブジェクトのタイルサイズを設定します。
     *
     * @param mb    新しい行のタイルサイズ (1より大きい必要があります)
     * @param nb    新しい列のタイルサイズ (1より大きい必要があります)
     */
    void set_ts(uint32_t mb, uint32_t nb);

    /**
     * @brief 行列の要素にランダムな数値を代入するメソッド。
     */
    void gen_rnd_elm();
    /**
     * @brief (ti,tj)要素へのポインタを返す。
     */
    T *elm(const uint32_t &ti, const uint32_t &tj) const;
    /**
     * @brief (ti,tj,i,j)要素へのポインタを返す。
     */
    T *elm(const uint32_t &ti, const uint32_t &tj,
           const uint32_t &i, const uint32_t &j) const;
    /**
     * @brief ファイル名を指定して行列の内容を出力するメソッド。
     * @param fname 出力先のファイル名。
     */
    void file_out(const char *fname);
    /**
     * @brief 行列の全要素をゼロに設定するメソッド。
     */
    void zero();
    /**
     * @brief 配列の各要素を加えるメソッド。
     * @param M 追加される行列。
     * @return 要素が追加された新しい行列。
     */
    Matrix<T> addElements(const Matrix<T> &M) const;
    /**
     * @brief 配列の各要素を引くメソッド。
     * @param M 減算される行列。
     * @return 要素が減算された新しい行列。
     */
    Matrix<T> minusElements(const Matrix<T> &M) const;

    /**
     * @brief 代入演算子のオーバーロード。
     * @param M 代入するMatrixオブジェクト。
     * @return 自身への参照。
     */
    Matrix &operator=(const Matrix &M);
    /**
     * @brief 加算演算子のオーバーロード。
     * @param M 加算するMatrixオブジェクト。
     * @return 加算結果の新しいMatrixオブジェクト。
     */
    Matrix operator+(const Matrix &M) const;
    /**
     * @brief 減算演算子のオーバーロード。
     * @param M 減算するMatrixオブジェクト。
     * @return 減算結果の新しいMatrixオブジェクト。
     */
    Matrix operator-(const Matrix &M) const;
    /**
     * @brief 等価演算子のオーバーロード。
     * @param M 比較元のMatrixオブジェクト。
     * @return 同等であればtrue、そうでなければfalseを返す。
     */
    bool operator==(const Matrix &M);
    /**
     * @brief 添字演算子のオーバーロード。
     * @param i アクセスする要素の位置。
     * @return i番目の要素への参照。
     */
    T &operator[](const uint64_t &i) const;
    /**
     * @brief 添字演算子のオーバーロード。
     * @param i 行の位置。
     * @param j 列の位置。
     * @return (i, j)位置の要素への参照。
     */
    T &operator()(const uint32_t &i, const uint32_t &j) const;
    /**
     * @brief 出力ストリームへの挿入演算子のフレンド関数。
     * @param os 出力ストリーム。
     * @param ma 表示するMatrix。
     * @return 行列の要素を出力したストリーム。
     */
    friend std::ostream& operator<< (std::ostream &os, const Matrix<T> &ma){
        uint32_t m = ma.m();
        uint32_t n = ma.n();

        for (uint32_t i = 0; i < m; ++i) {
            for (uint32_t j = 0; j < n; ++j) {
                os << ma(i, j) << " ";
            }
            os << std::endl;
        }

        return os;
    }

private:
    /**
 * @brief タイル行列から配列への変換を行う関数です。
 *
 * 与えられたタイル行列の座標 (`ti`, `tj`) および、その内部の座標 (`i`, `j`) を
 * 使用して、1次元配列に対応するインデックスを計算します。タイル行列から配列への
 * 変換は、タイル行列の配置 (Column-Major/Row-Major) とタイルの配置
 * (Column-Major/Row-Major) の組み合わせに基づく順序付けを制御します。
 *
 * @param ti タイル行列の行インデックスです。
 * @param tj タイル行列の列インデックスです。
 * @param i  タイル内の行インデックスです。
 * @param j  タイル内の列インデックスです。
 * @return 計算された1次元配列のインデックスを返します。
 */
    [[nodiscard]] uint64_t convertTileToArray(const uint32_t &ti, const uint32_t &tj,
                                              const uint32_t &i, const uint32_t &j) const;
};
