#ifndef __qwe_lwarr_h__
#define __qwe_lwarr_h__

template <class T> class LWArr {
        public: 
                LWArr() : m_n(0), m_a(0), m_p(0) {}
                LWArr(int n) : m_n(0), m_a(0), m_p(0) { rsz(n); }
                ~LWArr() { if (m_p) delete[](m_p); }
                int n() const { return m_n; }
                T* p() const { return m_p; }
                T* p(int i) { if (i>=m_n) resize(i+1); return m_p + i; }
                T* forget() { T* pp = m_p; m_n=m_a=0; m_p=0; return pp; }
                void clear() { if (m_p) delete[](m_p); m_n=m_a=0; m_p=0; }
                T* copy() { T* pp = new T[m_n]; memcpy(pp,m_p,sizeof(T)*m_n); return pp; }
                T& operator[] (int i) { return m_p[i]; }
                const T& operator[] (int i) const { return m_p[i]; }
                void add(const T& x) {
                        if (m_n==m_a) rsz(m_a ? 2*m_a : 4);
                        m_p[m_n++] = x;
                }
                T* add() {
                        if (m_n==m_a) rsz(m_a ? 2*m_a : 4);
                        return m_p + m_n++ ;
                }
                void ins_n(int i, int k) {
                        if (i>=m_n) { resize(i+k); return; }
                        resize(m_n+k);
                        memmove(m_p+i+k, m_p+i, sizeof(T) * (m_n-i-k) );
                }
                void ins(int i, const T& x) { ins_n(i,1); m_p[i]=x; }
                void cut(int i, int k=1) {
                        if (i>=m_n) return;
                        m_n-=k; if (i>=m_n) return;
                        memmove(m_p+i, m_p+i+k, sizeof(T) * (m_n-i) );
                }
                void add(const T& x, const T& y) { add(x); add(y); }
                void del(int i) { m_p[i] = m_p[--m_n]; }
                void resize(int n, bool rly=false) {
                        m_n = n;
                        if (rly) { rsz(n); return; }
                        if (n <= m_a) return;
                        for (n=(m_a?m_a*2:4); n<m_n; n*=2);
                        rsz(n);
                }
                void resize2(int n, const T& x, bool rly=false) {
                        int n0 = m_n;
                        resize(n,rly);
                        for (int i=n0; i<n; i++) m_p[i] = x;
                }
        protected:
                int m_n;
                int m_a;
                T * m_p;
                void rsz(int n) {
                        T* p = new T[n];
                        int n2 = n<m_n?n:m_n;
                        if (m_p) {
                                memcpy (p, m_p, n2*sizeof(T));
                                delete[] (m_p);
                        }
                        m_a = n;
                        m_p = p;
                }
};

#endif // __qwe_lwarr_h__
