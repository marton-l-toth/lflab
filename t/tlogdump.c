#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>

/* tlogdump (debug tool for lflab)
   dumps tlog in a more-or-less human-readable format (time:event pairs)
   time: duration (time until next event) in microseconds unless "m" says millisec.
   event: (s) - select (with wait)
   	  (S) - select (without wait)
   	  888 - calculate main mixer
          (p) - play (send to audio dev)
	  +-  - clock adjust
	  j88 - job slice (the longest one is displayed separately)
	  (J) - no time for job exec
   every 4th cycle the data is send to GUI process (for CPU meter display)
   (this has no separate event, you can see that every 4th "clock adjust" takes a little longer
*/

void dur4(char * to, int t) {
	if (t<=9999) return (void)sprintf(to, "%4d", t);
	if (t>262141) return (void)memcpy(to, "FRVR-WTF"+4*(t&1), 4);
	if (t>99999) return (void)sprintf(to, "%dm", t/1000);
	t/=100; sprintf(to, "%dm%c", t/10, 48+(t%10));
}

int ev4(char *to, int k) {
	if (k<128) to[0]=40, to[2]=41, to[3]=32, to[1] = k?k:63;
	else if (k<2048) return sprintf(to, "%+d", k-1090), k-1090;
	else if (k<4096) (k&=2047)<1000 ? sprintf(to, "j%d", k) 
				        : (k/=100, sprintf(to, "j%ck%c",48+k/10,48+k%10));
	else k&=4095, sprintf(to, "%d"+(k>999), k);
	return 0;
}

#define MEGA 1000000

int  ent10(char *q, unsigned int x) { return dur4(q,x&262143), q[4]=q[9]=58, ev4(q+5, x>>18); }
void top10(char *q, long long t) { sprintf(q, "%02d.%06d", (t/MEGA)%60, t%MEGA); }

int main(int ac, char ** av) {
	long long t0 = 0LL, t1 = 0LL;
	int nrow = 30, ncol = 10, eof = 0, nblk = 0, atot2 = 0;
	if (ac==4) nrow=atoi(av[1]), ncol = atoi(av[2]), av += 2;
	else if (ac!=2) return fprintf(stderr,"usage: %s [nrow ncol] <binfile>\n", *av), 1;
	int fd = open(av[1], O_RDONLY); if (fd<0) return  perror(av[1]), 1;
	int blks = nrow * ncol, ibuf[blks], wid = 10*ncol, osiz0 = (nrow+1)*wid;
	char obuf[osiz0 + 80]; memset(obuf, 32, osiz0);
	while(1) {
		int i,j,k,nr = read(fd, ibuf, 4*blks), eof = 4*blks - nr;
		t1 = t0;
		if (nr<=0) return nr ? (perror("read"),1) : 0;
		if (nr&3) fprintf(stderr, "ignoring %c extra bytes\n", nr&3);
		int nrow2 = ((nr>>=2) + ncol - 1) / ncol, osiz2 = (nrow2+1)*wid;
		int adjm = 0, adjM = 0, adjtot = 0;
		for (i=0,j=0,k=0; k<nr; k++, i += (++j>=nrow2 && !(j=0))) {
			if (!j) { top10(obuf+10*i, t0); obuf[10*i+9] = (i==ncol-1)?10:'|'; }
			int a = ent10(obuf + 10*i + wid*(j+1), ibuf[k]); t0 += ibuf[k]&262143;
			if (a) { adjtot+=a; if (a<adjm) adjm=a; else if (a>adjM) adjM=a; }
		}
		for (i=0; i<osiz2; i++) if (!obuf[i]) obuf[i] = 32;
		for (i=9; i<osiz2; i+=10) obuf[i] = 32; 
		for (i=wid-1; i<osiz2; i+=wid) obuf[i] = 10;
		int sec = t0/MEGA, nx;
		nx=sprintf(obuf+osiz2, "blk%04d t%02d:%02d.%06d len:%d.%06d +-min,max,tot,atot: %d,%d,%d,%d\n",
		    nblk++,sec/60,sec%60,t0%MEGA,(t0-t1)/MEGA,(t0-t1)%MEGA,adjm,adjM,adjtot,atot2+=adjtot);
		write(1, obuf, osiz2+nx);
		if (eof) return 0;
	}}
