#ifndef __qwe_aobuf_h__
#define __qwe_aobuf_h__

class AOBuf {
        public: 
		AOBuf() : m_ec(0) {}
		virtual ~AOBuf() {}
                virtual int vpf(const char * fmt, va_list ap) = 0;
                virtual int sn(const char * s, int n) = 0;
                virtual void flush() {}
		int pf(const char * fmt, ...) { 
			va_list ap; va_start(ap, fmt);
			int k = vpf(fmt, ap); va_end(ap); return k; 
		}
		int err() { int r = m_ec; m_ec = 0; return r; }
	protected:
		int m_ec;
};

class FOBuf : public AOBuf {
        public: 
                FOBuf(FILE *f) : m_f(f) {}
		virtual ~FOBuf() { if (m_f) fclose(m_f);  }
                FILE * f() const { return m_f; }
                virtual int vpf(const char * fmt, va_list ap) {
			if (m_ec) return m_ec;
			int k = vfprintf(m_f, fmt, ap); return (k<0) ? (m_ec = EEE_ERRNO) : k;
		}
                virtual int sn(const char * s, int n) { 
			if (m_ec) return m_ec;
			int k = fwrite(s, 1, n, m_f); return (k<0) ? (m_ec = EEE_ERRNO) : k;
		}
                virtual void flush() { fflush(m_f); }
        protected:
                FILE * m_f;
};

static int xprintf(AOBuf *p, const char * fmt, ...) {
	va_list ap; va_start(ap, fmt);
	int k = p->vpf(fmt, ap); va_end(ap); return k; 
}

#endif // __qwe_aobuf_h__
