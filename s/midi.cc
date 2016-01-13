#include <math.h>
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

#include "uc1.h"
#include "glob.h"

int rawmidi_desc(char *to, int id, int maxlen); //asnd.cc

class AMidi {
	public:
		virtual void cmd0(int x) = 0;
		virtual void cmd1(int x, int y) = 0;
		virtual void cmd2(int x, int y, int z) = 0;
		virtual void longcmd(const unsigned char *p, int n) = 0;
		const char * nm() { return m_nm; }
		int id() const { return m_id; }
		int ini(int id);
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

int AMidi::ini(int id) {
	if (id>=0) m_id = id; else id = m_id;
	rawmidi_desc(m_nm, id, 64);
	devname[14] = 48+(id>>4); devname[16] = 48+(id&15);
	int fd = open(devname, O_RDONLY); 
	if (fd<0) log("midi: %s(%s): %s", devname, m_nm, strerror(errno));
	else 	  log("midi: %s(%s): opened(%d)", devname , m_nm, fd);
	return fd;
}

void LogMidi::longcmd(const unsigned char * p, int n) {
	log_n("midi-c+(%d,%s) : [", m_id, m_nm);
	for (int i=0; i<n; i++) log_n(" %0x", p[i]); log(" ]"); }

unsigned int midi_bv;
int midi_fd[32];
static AMidi * midi_p[32];

void midi_input(int i) {
	unsigned char buf[256], *p = buf;
	int r = read(midi_fd[i], buf, 256); if (r<0) perror(midi_p[i]->nm());
	AMidi * q = midi_p[i];
	while (r) {
		int x = *p;
		if (x<0xf0) {
			if ((x&0xe0)==0xc0) { if (r<2) return log("midi-frag1: %s(0x%x)", q->nm(), q->id());
					      q->cmd1(x, p[1]); p+=2; r -= 2; }
			else if (r<3) return log("midi-frag%c: %s(0x%x)", r+48, q->nm(), q->id());
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
	for (int i=0; i<n; i++) if ((midi_fd[i] = (midi_p[i] = new LogMidi())->ini(buf[i])) < 0) 
		delete(midi_p[i]); else midi_bv |= 1u<<i;
}

