#ifndef __VISUALIZE_H
#define __VISUALIZE_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Imlib2.h>

template <typename T, int _Flags>
void visualize(Eigen::SparseMatrix<T, _Flags> &matrix, const char* filename) {
    int my = matrix.rows();
    int mx = matrix.cols();
    if (mx * my > 1000000) {
        std::cerr <<
            "visualize(): warning: matrix too large, not writing "
            "image `" << filename << "'\n";
        return;
    }
    unsigned char *pixmap = (unsigned char*)malloc(4 * mx * my);
    memset(pixmap, 255,  4 * mx * my);
    for (int k = 0; k < matrix.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(matrix,k); it; ++it) {
            int i = it.col() * mx + it.row();
            pixmap[4 * i    ] = 0;
            pixmap[4 * i + 1] = 0;
            pixmap[4 * i + 2] = 0;
        }
    }
    Imlib_Image img = imlib_create_image_using_data(mx, my, (DATA32 *) pixmap);
    imlib_context_set_image(img);
    imlib_save_image(filename);
}


// vim: ts=4 sw=4 expandtab
#endif /* __VISUALIZE_H */
