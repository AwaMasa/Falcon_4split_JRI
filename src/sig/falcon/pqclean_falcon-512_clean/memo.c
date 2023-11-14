// ２分割→４分割...
void
PQCLEAN_FALCON512_CLEAN_poly_split_fft(
    fpr *f0, fpr *f1, fpr *f2, fpr *f3,
    const fpr *f, unsigned logn) {
    /*
     * The FFT representation we use is in bit-reversed order
     * (element i contains f(w^(rev(i))), where rev() is the
     * bit-reversal function over the ring degree. This changes
     * indexes with regards to the Falcon specification.
     */
    size_t n, hn, qn, on, u;

    n = (size_t)1 << logn;
    hn = n >> 1;
    qn = hn >> 1;
    on = qn >> 1;

    /*
     * We process complex values by pairs. For logn = 2, there is only
     * one complex value (the other one is the implicit conjugate),
     * so we add the two lines below because the loop will be
     * skipped.
     */
    f0[0] = f[0];
    f1[0] = f[qn];
    f2[0] = f[hn];
    f3[0] = f[qn*3];
    
    // for (u = 0; u <= hn; u+=2) {
    for (u = 0; u < on; u ++) {
        fpr a_re, a_im, b_re, b_im;
        fpr c_re, c_im, d_re, d_im;

        fpr ab_add_re, ab_add_im, ab_sub_re, ab_sub_im;
        fpr cd_add_re, cd_add_im, cd_sub_re, cd_sub_im;
        fpr t_re, t_im;

        a_re = f[(u << 2) + 0];
        a_im = f[(u << 2) + 0 + hn];
        b_re = f[(u << 2) + 1];
        b_im = f[(u << 2) + 1 + hn];
        c_re = f[(u << 2) + 2];
        c_im = f[(u << 2) + 2 + hn];
        d_re = f[(u << 2) + 3];
        d_im = f[(u << 2) + 3 + hn];


        // gm[16] → (((u≪1)+hn)≪1)+0
        // gm[17] → (((u≪1)+hn)≪1)+1
        // gm[18] → (((u≪1)+hn)≪1)+2
        // gm[19] → (((u≪1)+hn)≪1)+3

        // gm[8 ]→ ((u+qn)≪1)+0
        // gm[9 ]→ ((u+qn)≪1)+1

        FPC_ADD(ab_add_re, ab_add_im, a_re, a_im, b_re, b_im);

        FPC_SUB(ab_sub_re, ab_sub_im, a_re, a_im, b_re, b_im);
        // FPC_MUL(ab_sub_re, ab_sub_im, ab_sub_re, ab_sub_im,
        //     fpr_gm_tab[(((u << 1) + hn) << 1) + 0],
        //     fpr_neg(fpr_gm_tab[(((u << 1) + hn) << 1) + 1]));

        FPC_ADD(cd_add_re, cd_add_im, c_re, c_im, d_re, d_im);

        FPC_SUB(cd_sub_re, cd_sub_im, c_re, c_re, d_re, d_im);
        // FPC_MUL(cd_sub_re, cd_sub_im, cd_sub_re, cd_sub_im,
        //     fpr_gm_tab[(((u << 1) + hn) << 1) + 2],
        //     fpr_neg(fpr_gm_tab[(((u << 1) + hn) << 1) + 3]));


        // f0 
        FPC_ADD(t_re, t_im, ab_add_re, ab_add_im, cd_add_re, cd_add_im);
        f0[u] = fpr_quarter(t_re);
        f0[u + on] = fpr_quarter(t_im);

        // f1 
        FPC_SUB(t_re, t_im, ab_add_re, ab_add_im, cd_add_re, cd_add_im);
        FPC_MUL(t_re, t_im, t_re, t_im,
            fpr_gm_tab[((u + qn) << 1) + 0],
            fpr_neg(fpr_gm_tab[((u + qn) << 1) + 1]));
        f1[u] = fpr_quarter(t_re);
        f1[u + on] = fpr_quarter(t_im);

        // a_re <- (F_0−F_8 )gm[16]−(F_1−F_9 )gm[17]
        // a_im <- (F_0−F_8 )gm[17]+(F_1−F_9 )gm[16]　
        FPC_MUL(a_re, a_im, ab_sub_re, ab_sub_im,
            fpr_gm_tab[(((u << 1) + hn) << 1) + 0],
            fpr_neg(fpr_gm_tab[(((u << 1) + hn) << 1) + 1]));

        // b_re <- (F_4−F_12 )gm[18]−(F_5−F_13 )gm[19]
        // b_im <- (F_4−F_12 )gm[19]+(F_5−F_13 )gm[18]
        FPC_MUL(b_re, b_im, cd_sub_re, cd_sub_im,
            fpr_gm_tab[(((u << 1) + hn) << 1) + 2],
            fpr_neg(fpr_gm_tab[(((u << 1) + hn) << 1) + 3]));

        // f2
        FPC_ADD(t_re, t_im, a_re, a_im, b_re, b_im);
        f2[u] = fpr_quarter(t_re);
        f2[u + on] = fpr_quarter(t_im);

        // c_re <- ((F_0−F_8 )gm[16]−(F_1−F_9 )gm[17]) - ((F_4−F_12 )gm[18]−(F_5−F_13 )gm[19])
        // c_im <- ((F_0−F_8 )gm[17]+(F_1−F_9 )gm[16]) - ((F_4−F_12 )gm[19]+(F_5−F_13 )gm[18])
        FPC_SUB(c_re, c_im, a_re, a_im, b_re, b_im);
        

        // f3
        FPC_MUL(t_re, t_im, c_re, c_im,
            fpr_gm_tab[((u + qn) << 1) + 0],
            fpr_neg(fpr_gm_tab[((u + qn) << 1) + 1]));
        f3[u] = fpr_quarter(t_re);
        f3[u + on] = fpr_quarter(t_im);


    }
}