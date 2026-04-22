#ifndef SRC_HPP
#define SRC_HPP

#include <vector>
#include "fraction.hpp"

class matrix {
private:
    int m, n;
    fraction **data;
    void alloc(int rows, int cols){ m=rows; n=cols; if(m==0||n==0){ data=nullptr; return; } data=new fraction*[m]; for(int i=0;i<m;++i) data[i]=new fraction[n]; }
    void free_data(){ if(!data) return; for(int i=0;i<m;++i) delete [] data[i]; delete [] data; data=nullptr; m=n=0; }
public:
    matrix(): m(0), n(0), data(nullptr) {}
    matrix(int m_, int n_): m(0), n(0), data(nullptr){ if(m_<0||n_<0) throw matrix_error(); alloc(m_, n_); }
    matrix(const matrix& o): m(0), n(0), data(nullptr){ alloc(o.m,o.n); for(int i=0;i<m;++i) for(int j=0;j<n;++j) data[i][j]=o.data[i][j]; }
    matrix(matrix&& o) noexcept: m(o.m), n(o.n), data(o.data){ o.m=o.n=0; o.data=nullptr; }
    ~matrix(){ free_data(); }
    matrix& operator=(const matrix& o){ if(this==&o) return *this; free_data(); alloc(o.m,o.n); for(int i=0;i<m;++i) for(int j=0;j<n;++j) data[i][j]=o.data[i][j]; return *this; }
    fraction& operator()(int i, int j){ if(i<=0||i>m||j<0||j>=n||!data) throw matrix_error(); return data[i-1][j]; }
    friend matrix operator*(const matrix& a, const matrix& b){ if(!a.data||!b.data||a.n!=b.m) throw matrix_error(); matrix r(a.m,b.n); for(int i=0;i<a.m;++i){ for(int k=0;k<a.n;++k){ fraction aik=a.data[i][k]; if(aik==fraction(0)) continue; for(int j=0;j<b.n;++j) r.data[i][j]=r.data[i][j]+aik*b.data[k][j]; } } return r; }
    matrix transposition(){ if(!data||m==0||n==0) throw matrix_error(); matrix t(n,m); for(int i=0;i<m;++i) for(int j=0;j<n;++j) t.data[j][i]=data[i][j]; return t; }
    fraction determination(){ if(!data||m==0||n==0||m!=n) throw matrix_error(); matrix tmp(*this); fraction det(1); int sign=1; for(int c=0;c<n;++c){ int piv=-1; for(int r=c;r<m;++r) if(!(tmp.data[r][c]==fraction(0))){ piv=r; break; } if(piv==-1) return fraction(0); if(piv!=c){ auto row=tmp.data[piv]; tmp.data[piv]=tmp.data[c]; tmp.data[c]=row; sign=-sign; } fraction pv=tmp.data[c][c]; det=det*pv; for(int r=c+1;r<m;++r){ if(tmp.data[r][c]==fraction(0)) continue; fraction f=tmp.data[r][c]/pv; for(int j=c;j<n;++j) tmp.data[r][j]=tmp.data[r][j]-f*tmp.data[c][j]; } } if(sign==-1) det=fraction(0)-det; return det; }
    static std::vector<fraction> solve_linear(matrix A, const std::vector<fraction>& b){ if(!A.data||A.m==0||A.n==0||(int)b.size()!=A.m) throw matrix_error(); int rows=A.m, cols=A.n; std::vector<std::vector<fraction>> aug(rows, std::vector<fraction>(cols+1)); for(int i=0;i<rows;++i){ for(int j=0;j<cols;++j) aug[i][j]=A.data[i][j]; aug[i][cols]=b[i]; } int r=0; for(int c=0;c<cols && r<rows;++c){ int piv=-1; for(int i=r;i<rows;++i) if(!(aug[i][c]==fraction(0))){ piv=i; break; } if(piv==-1) continue; if(piv!=r) std::swap(aug[piv], aug[r]); fraction pv=aug[r][c]; for(int j=c;j<=cols;++j) aug[r][j]=aug[r][j]/pv; for(int i=0;i<rows;++i){ if(i==r) continue; if(aug[i][c]==fraction(0)) continue; fraction f=aug[i][c]; for(int j=c;j<=cols;++j) aug[i][j]=aug[i][j]-f*aug[r][j]; } ++r; } std::vector<fraction> x(cols, fraction(0)); for(int i=0;i<rows;++i){ int lead=-1; for(int j=0;j<cols;++j) if(!(aug[i][j]==fraction(0))){ lead=j; break; } if(lead!=-1) x[lead]=aug[i][cols]; } return x; }
};

class resistive_network{
private:
    int interface_size, connection_size;
    matrix adjacency, conduction;
public:
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[])
        : interface_size(interface_size_), connection_size(connection_size_), adjacency(interface_size_, interface_size_), conduction(interface_size_, interface_size_){
        for(int e=0;e<connection_size_;++e){ int u=from[e]-1, v=to[e]-1; fraction g=fraction(1)/resistance[e]; adjacency(u+1,v)=adjacency(u+1,v)+g; adjacency(v+1,u)=adjacency(v+1,u)+g; conduction(u+1,u)=conduction(u+1,u)+g; conduction(v+1,v)=conduction(v+1,v)+g; conduction(u+1,v)=conduction(u+1,v)-g; conduction(v+1,u)=conduction(v+1,u)-g; }
    }
    fraction get_equivalent_resistance(int interface_id1, int interface_id2){ int n=interface_size; std::vector<fraction> b(n, fraction(0)); b[interface_id1-1]=b[interface_id1-1]+fraction(1); b[interface_id2-1]=b[interface_id2-1]-fraction(1); int dim=n-1; matrix C(dim,dim); std::vector<fraction> br(dim); for(int i=0;i<dim;++i){ for(int j=0;j<dim;++j) C(i+1,j)=conduction(i+1,j); br[i]=b[i]; } std::vector<fraction> u=matrix::solve_linear(C, br); fraction v1=(interface_id1==n)?fraction(0):u[interface_id1-1]; fraction v2=(interface_id2==n)?fraction(0):u[interface_id2-1]; return v1 - v2; }
    fraction get_voltage(int id, fraction current[]){ int n=interface_size; int dim=n-1; matrix C(dim,dim); std::vector<fraction> br(dim); for(int i=0;i<dim;++i){ for(int j=0;j<dim;++j) C(i+1,j)=conduction(i+1,j); br[i]=current[i]; } std::vector<fraction> u=matrix::solve_linear(C, br); if(id==n) return fraction(0); return u[id-1]; }
    fraction get_power(fraction voltage[]){ int n=interface_size; std::vector<fraction> y(n, fraction(0)); for(int i=0;i<n;++i){ fraction sum(0); for(int j=0;j<n;++j) sum=sum+conduction(i+1,j)*voltage[j]; y[i]=sum; } fraction P(0); for(int i=0;i<n;++i) P=P+voltage[i]*y[i]; return P; }
};

#endif
