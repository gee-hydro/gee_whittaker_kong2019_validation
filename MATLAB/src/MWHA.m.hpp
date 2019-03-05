// Automatically translated using m2cpp 2.0 on 2019-03-03 12:32:17

#ifndef MWHA_M_HPP
#define MWHA_M_HPP

#include "mconvert.h"
#include <armadillo>
using namespace arma ;

void MWHA(vec x, vec xi, int nf, int dm, int HiLo, double thr, double fet, double high, double low, vec yi, vec& y_1, vec& y_2, vec& y_3, vec& y_or)
{
  TYPE diff11, diff12, diff13, diff14, diff21, diff22, diff31, diff32, diff_yi, lo1, lo_py2, mat_mean1, mat_mean2, mat_mean3, p1, p2, p3, p4, p5, p6, p_y_or1, p_y_or2, p_y_or3, p_y_or4, p_y_or5, p_yi1, p_yi2, p_yi3, p_yi4, p_yi5, p_yi6, rankVec, temp, w_new, w_or ;
  int dm2, dmi, i_num, j, m, noutmax, wi, wj ;
  mat A, B, V, YY, d1, d2, di, p, px, v, w ;
  rowvec p_num, p_yi, r ;
  uvec lo, lo2, lo3, lo_diff1, lo_p1, lo_p2, lo_p3, lo_p4, lo_p5, lo_p6, lo_py1, num_0 ;
  uword nnodes, nout, npoints ;
  vec Y_mw, ang, cs1, cs2, diffVec, max_y_err, sn1, sn2, x1, y_err, yi1, yi2 ;
  nnodes = m2cpp::length(xi) ;
  num_0 = find(yi, 1) + 1 ;
  if (m2cpp::isempty(num_0))
  {
    y_or = arma::zeros<vec>(nnodes) ;
  }
  else
  {
    p_yi = arma::ones<rowvec>(nnodes) ;
    diff_yi = diff(yi) ;
    dm2 = 10 ;
    if (HiLo==1)
    {
      lo_py1 = find(diff_yi<=-thr) + 1 ;
      lo_py2 = diff_yi>=thr ;
      p_yi(lo_py1).fill(0) ;
      p_yi(lo_py2-1).fill(0) ;
    }
    else if (HiLo==-1)
    {
      lo_py1 = find(diff_yi>=thr) + 1 ;
      lo_py2 = diff_yi<=-thr ;
      p_yi(lo_py1).fill(0) ;
      p_yi(lo_py2-1).fill(0) ;
    }
    lo = find(p_yi==0) + 1 ;
    x1 = x(lo-1) ;
    npoints = m2cpp::length(x1) ;
    V = arma::zeros<mat>(2*nf+1, npoints) ;
    w = arma::zeros<mat>(npoints, nnodes) ;
    p = arma::ones<mat>(2*nf+1, nnodes) ;
    px = arma::ones<mat>(2*nf+1, npoints) ;
    noutmax = nnodes-2*nf-2 ;
    ang = arma::trans((2*datum::pi*(m2cpp::fspan(1, 1, nf))/nnodes)) ;
    cs1 = arma::cos(ang*xi) ;
    sn1 = arma::sin(ang*xi) ;
    cs2 = arma::cos(ang*x1) ;
    sn2 = arma::sin(ang*x1) ;
    for (m=1; m<=nf; m++)
    {
      p.row(2*m-1) = arma::strans(cs1(m, m2cpp::span<uvec>(0, cs1.n_cols-1))) ;
      p.row(2*m) = arma::strans(sn1(m, m2cpp::span<uvec>(0, sn1.n_cols-1))) ;
      px.row(2*m-1) = arma::strans(cs2(m, m2cpp::span<uvec>(0, cs2.n_cols-1))) ;
      px.row(2*m) = arma::strans(sn2(m, m2cpp::span<uvec>(0, sn2.n_cols-1))) ;
    }
    d1 = arma::ones<mat>(nnodes, npoints)*diagmat(x1) ;
    d2 = arma::ones<mat>(npoints, nnodes)*diagmat(xi) ;
    di = arma::strans(d1)-d2 ;
    p_num = arma::ones<rowvec>(nnodes) ;
    p_num(yi<low || yi>high-1).fill(0) ;
    p_num(yi==0-1).fill(0) ;
    for (j=1; j<=npoints; j++)
    {
      dmi = dm2 ;
      nout = nnodes ;
      while (nout>noutmax && dmi<5*dm2)
      {
        r = abs(di.row(j-1))/dmi ;
        lo1 = r>1.0 ;
        lo2 = find(r<=0.5) + 1 ;
        lo3 = find(r>0.5 && r<=1) + 1 ;
        wj(lo1) = 0 ;
        wj(lo2) = int(arma::as_scalar(2/3.0-4*arma::square(arma::strans(r(lo2-1)))+4*arma::pow(arma::strans(r(lo2-1)), 3))) ;
        wj(lo3) = int(arma::as_scalar(4/3.0-4*arma::strans(r(lo3-1))+4*arma::square(arma::strans(r(lo3-1)))-4*arma::pow(arma::strans(r(lo3-1)), 3)/3.0)) ;
        w.row(j-1) = wj*p_num*p_yi ;
        nout = arma::sum(w.row(j-1)==0) ;
        dmi = dmi+1 ;
      }
      if ((nout>noutmax))
      {
        return ;
      }
      B = p*diagmat(w.row(j-1)) ;
      A = B*arma::trans(p) ;
      v = arma::solve(A, B, solve_opts::fast)*yi ;
      V.col(j-1) = v ;
    }
    YY = arma::trans(V)*px ;
    Y_mw = diagvec(YY) ;
    yi(lo-1) = Y_mw ;
    npoints = m2cpp::length(x) ;
    V = arma::zeros<mat>(2*nf+1, npoints) ;
    w = arma::zeros<mat>(npoints, nnodes) ;
    p = arma::ones<mat>(2*nf+1, nnodes) ;
    px = arma::ones<mat>(2*nf+1, npoints) ;
    noutmax = nnodes-2*nf-2 ;
    ang = arma::trans((2*datum::pi*(m2cpp::fspan(1, 1, nf))/nnodes)) ;
    cs1 = arma::cos(ang*xi) ;
    sn1 = arma::sin(ang*xi) ;
    cs2 = arma::cos(ang*x) ;
    sn2 = arma::sin(ang*x) ;
    for (m=1; m<=nf; m++)
    {
      p.row(2*m-1) = arma::strans(cs1(m, m2cpp::span<uvec>(0, cs1.n_cols-1))) ;
      p.row(2*m) = arma::strans(sn1(m, m2cpp::span<uvec>(0, sn1.n_cols-1))) ;
      px.row(2*m-1) = arma::strans(cs2(m, m2cpp::span<uvec>(0, cs2.n_cols-1))) ;
      px.row(2*m) = arma::strans(sn2(m, m2cpp::span<uvec>(0, sn2.n_cols-1))) ;
    }
    d1 = arma::ones<mat>(nnodes, npoints)*diagmat(x) ;
    d2 = arma::ones<mat>(npoints, nnodes)*diagmat(xi) ;
    di = arma::strans(d1)-d2 ;
    p_num = arma::ones<rowvec>(nnodes) ;
    p_num(yi<low || yi>high-1).fill(0) ;
    p_num(yi==0-1).fill(0) ;
    y_or = yi ;
    max_y_err = fet+1 ;
    i_num = 0 ;
    while (max_y_err>fet && i_num<1000)
    {
      i_num = i_num+1 ;
      for (j=1; j<=npoints; j++)
      {
        dmi = dm ;
        nout = nnodes ;
        while (nout>noutmax && dmi<5*dm)
        {
          r = abs(di.row(j-1))/dmi ;
          lo1 = r>1.0 ;
          lo2 = find(r<=0.5) + 1 ;
          lo3 = find(r>0.5 && r<=1) + 1 ;
          wi(lo1) = 0 ;
          wi(lo2) = int(arma::as_scalar(2/3.0-4*arma::square(arma::strans(r(lo2-1)))+4*arma::pow(arma::strans(r(lo2-1)), 3))) ;
          wi(lo3) = int(arma::as_scalar(4/3.0-4*arma::strans(r(lo3-1))+4*arma::square(arma::strans(r(lo3-1)))-4*arma::pow(arma::strans(r(lo3-1)), 3)/3.0)) ;
          w.row(j-1) = wi*p_num ;
          nout = arma::sum(w.row(j-1)==0) ;
          dmi = dmi+1 ;
        }
        if ((nout>noutmax))
        {
          return ;
        }
        B = p*diagmat(w.row(j-1)) ;
        A = B*arma::trans(p) ;
        v = arma::solve(A, B, solve_opts::fast)*y_or ;
        V.col(j-1) = v ;
      }
      YY = arma::trans(V)*px ;
      Y_mw = diagvec(YY) ;
      p_num = arma::ones<rowvec>(nnodes) ;
      if (i_num==1)
      {
        y_1 = Y_mw ;
      }
      if (HiLo==1)
      {
        y_err = y_or-Y_mw ;
        diffVec = arma::trans(p_num)%y_err ;
        lo_diff1 = find(diffVec<0) + 1 ;
        y_or(lo_diff1-1) = Y_mw(lo_diff1-1) ;
        [~, rankVec] = sort(diffVec, "ascend") ;
        max_y_err = y_err(rankVec(nnodes)-1) ;
      }
      else
      {
        max_y_err = fet-1 ;
        y_or = Y_mw ;
      }
      if (i_num==1 || i_num==2)
      {
        y_2 = y_or ;
      }
      p_num = arma::ones<rowvec>(nnodes) ;
    }
    y_3 = y_or ;
    yi1 = y_1 ;
    yi2 = yi ;
    mat_mean1 = arma::ones<vec>(nnodes)*mean(yi) ;
    mat_mean2 = arma::ones<vec>(nnodes)*mean(yi(yi>mean(yi)-1)) ;
    mat_mean3 = arma::ones<vec>(nnodes)*mean(yi(yi<mean(yi)-1)) ;
    diff11 = y_or-mat_mean1 ;
    diff12 = yi-mat_mean1 ;
    diff13 = mat_mean2-y_or ;
    diff14 = mat_mean2-yi ;
    p_y_or2 = diff11>=0 ;
    p_yi2 = diff12>=0 ;
    p_y_or3 = diff13>=0 ;
    p_yi3 = diff14>=0 ;
    p2 = p_y_or2%p_yi2%p_y_or3%p_yi3 ;
    lo_p2 = find(p2==1) + 1 ;
    w_new = (diff11-diff12)/diff11 ;
    w_or = diff12/diff11 ;
    y_or(lo_p2-1) = w_new(lo_p2)%y_or(lo_p2-1)+w_or(lo_p2)%yi2(lo_p2-1) ;
    diff21 = y_or-mat_mean2 ;
    diff22 = yi-mat_mean2 ;
    p_y_or1 = diff21>=0 ;
    p_yi1 = diff22>=0 ;
    p1 = p_y_or1%p_yi1 ;
    lo_p1 = find(p1==1) + 1 ;
    w_new = (diff21-diff22)/diff21 ;
    w_or = diff22/diff21 ;
    y_or(lo_p1-1) = w_new(lo_p1)%y_or(lo_p1-1)+w_or(lo_p1)%yi2(lo_p1-1) ;
    diff31 = y_or-mat_mean3 ;
    diff32 = yi-mat_mean3 ;
    p_y_or4 = diff11<=0 ;
    p_yi4 = diff12<=0 ;
    p_y_or5 = diff31>=0 ;
    p_yi5 = diff32>=0 ;
    p_yi6 = diff32<=0 ;
    p3 = p_y_or4%p_yi4%p_y_or5%p_yi5 ;
    lo_p3 = find(p3==1) + 1 ;
    w_new = (diff31-diff32)/diff31 ;
    w_or = diff32/diff31 ;
    y_or(lo_p3-1) = w_new(lo_p3)%y_or(lo_p3-1)+w_or(lo_p3)%yi2(lo_p3-1) ;
    p4 = p_y_or1%p_yi3%p_yi2 ;
    lo_p4 = find(p4==1) + 1 ;
    w_new = diff21/(diff21+diff14) ;
    w_or = diff14/(diff21+diff14) ;
    temp = w_new ;
    w_new = max(w_new, w_or) ;
    w_or = min(temp, w_or) ;
    y_or(lo_p4-1) = w_new(lo_p4)%y_or(lo_p4-1)+w_or(lo_p4)%yi1(lo_p4-1) ;
    p5 = p_y_or2%p_yi4%p_yi5%p_y_or4 ;
    lo_p5 = find(p5==1) + 1 ;
    w_new = diff11/(diff11+abs(diff12)) ;
    w_or = abs(diff12)/(diff11+abs(diff12)) ;
    temp = w_new ;
    w_new = max(w_new, w_or) ;
    w_or = min(temp, w_or) ;
    y_or(lo_p5-1) = w_new(lo_p5)%y_or(lo_p5-1)+w_or(lo_p5)%yi1(lo_p5-1) ;
    p6 = p_y_or5%p_yi6%p_y_or4 ;
    lo_p6 = find(p6==1) + 1 ;
    w_new = diff31/(diff31+abs(diff32)) ;
    w_or = abs(diff32)/(diff31+abs(diff32)) ;
    temp = w_new ;
    w_new = max(w_new, w_or) ;
    w_or = min(temp, w_or) ;
    y_or(lo_p6-1) = w_new(lo_p6)%y_or(lo_p6-1)+w_or(lo_p6)%yi1(lo_p6-1) ;
  }
}
#endif