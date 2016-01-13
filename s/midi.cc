#include <math.h>
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

#include "uc1.h"
#include "glob.h"
#include "errtab.inc"

int rawmidi_desc(char *to, int id, int maxlen); //asnd.cc

class AMidi {
	public:
		virtual ~AMidi() {}
		virtual void cmd0(int x) = 0;
		virtual void cmd1(int x, int y) = 0;
		virtual void cmd2(int x, int y, int z) = 0;
		virtual void longcmd(const unsigned char *p, int n) = 0;
		const char * nm() { return m_nm; }
		int id() const { return m_id; }
		int ini(int ix, int id);
		void errmsg(int ec);
	protected:
		int m_id;
		char m_nm[64];
};

class LogMidi : public AMidi {
	public:
		virtual void cmd0(int x) { log("midi_c0(%d,%s) : %02x", m_id, m_nm, x); }
		virtual void cmd1(int x,int y) { log("midi_c1(%d,%s) : %02x %02x", m_id, m_nm,x,y);}
		virtual void cmd2(int x,int y,int z) { log("midi_c2(%d,%s) : %02x %02x %02x",
				m_id, m_nm, x, y, z); }
		virtual void longcmd(const unsigned char *p, int n);
};

static char devname[24];

int AMidi::ini(int ix, int id) {
	if (id>=0) m_id = id; else id = m_id;
	rawmidi_desc(m_nm, id, 64);
	devname[14] = 48+(id>>4); devname[16] = 48+(id&15);
	int fd = open(devname, (debug_flags&DFLG_RWMIDI) ? O_RDWR : O_RDONLY); 
	if (fd<0) log("midi%x: %s(%s): %s", ix, devname, m_nm, strerror(errno));
	else 	  log("midi%x: %s(%s): opened(%d)", ix, devname , m_nm, fd);
	return fd;
}

void AMidi::errmsg(int ec) { gui_errq_add(ec); log("midi (0x%x:%s): %s", m_id, m_nm, err_str(ec)); }

void LogMidi::longcmd(const unsigned char * p, int n) {
	log_n("midi-c+(%d,%s) : [", m_id, m_nm);
	for (int i=0; i<n; i++) log_n(" %0x", p[i]); log(" ]"); }

unsigned int midi_bv;
int midi_fd[32];
static AMidi * midi_p[32];

int midi_cmd(const char *s) {
	if (!s || !*s) return MDE_PARSE;
	unsigned char buf[64];
	int ec, ix = b32_to_i(*(s++)), n = 0;
	for (int i=0; i<64; i++) {
		while (*s==32) ++s;
		if (!*s) break;
		if (!s[1]) return MDE_PARSE;
		buf[n++] = hex2(s); s += 2;
	}
	if ((ec=write(midi_fd[ix], buf, n))<=0) return ec ? EEE_ERRNO : EEE_ZEROLEN;
	return 0;
}


void midi_input(int i) {
	unsigned char buf[256], *p = buf;
	AMidi * q = midi_p[i];
	int r = read(midi_fd[i], buf, 256); 
	if (r<0) return close(midi_fd[i]), midi_p[i]=0, midi_bv &= ~(1<<i), // TODO: reopen for 0 (?)
			q->errmsg(EEE_ERRNO), q->errmsg(MDE_FATAL), delete(q);
	while (r) {
		int x = *p;
		if (x<0xf0) {
			if (x<128) return q->errmsg(MDE_SEVEN);
			if ((x&0xe0)==0xc0) { if (r<2 || p[1]>127) return q->errmsg(MDE_FRAG12);
					      q->cmd1(x, p[1]); p+=2; r -= 2; }
			else if (r<3 || ((p[1]|p[2])&128) ) return q->errmsg(MDE_FRAG23-(r==1||p[2]>127));
			q->cmd2(x, p[1], p[2]); p+=3, r-=3;
			continue;
		}
		if (x==0xf0) {
			int j = 1; while (j<r && p[j]!=0xf7) j++;
			if (j==r) return log("midi-frag%d: %s(0x%x)", r, q->nm(), q->id());
			q->longcmd(p+1, j-1); p += j+1, r -= j+1;
			continue;
		}
		if (x>0xf3) q->cmd0(x), p++, r--;
		else (x==0xf2) ? (q->cmd2(x, p[1], p[2]), p+=3, r-=3) : (q->cmd1(x, p[1]), p+=2, r-=2);
	}
}

void midi_init() {
	memcpy(devname, "/dev/snd/midiCxDy", 18);
	unsigned char buf[32];
	int n = find_dev(buf, 1, 32);
	for (int i=0; i<n; i++) if ((midi_fd[i] = (midi_p[i] = new LogMidi())->ini(i,buf[i])) < 0) 
		delete(midi_p[i]); else midi_bv |= 1u<<i;
}

