#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <limits.h>
#include <sys/select.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/inotify.h>

#define FP0(S) 		fprintf(stderr, S "\n")
#define FP1(S,X) 	fprintf(stderr, S "\n", (X))
#define FP2(S,X,Y) 	fprintf(stderr, S "\n", (X),(Y))
#define FP3(S,X,Y,Z) 	fprintf(stderr, S "\n", (X),(Y),(Z))
#define FP4(S,X,Y,Z,W) 	fprintf(stderr, S "\n", (X),(Y),(Z),(W))

static const char * ino_dir;

void cmd(int ino, char *s, int l) {
	int k; s[l-(l&&s[l-1]==10)] = 0;
	switch(*s) {
		case 0  : FP0("bye"), rmdir(ino_dir), exit(0);
		default:  FP0("hmm???"); return;
	}}

#define INEV_SIZ (sizeof(struct inotify_event))

void evh(struct inotify_event *ev) {
	FP2("msk: 0x%x, nm=\"%s\"", ev->mask, ev->name);
}

void hee(int ino) {
	char buf[4096], *s = buf;   int r;
	if ((r = read(ino, buf, 4096))<=0) { perror("read/ino"); return; }
	while(1) {
		struct inotify_event *ev = (struct inotify_event*)s;
		evh(ev);
		int l = INEV_SIZ + ev->len; if ((r-=l)<=0) return; else s+=l;
	}}

int main(int ac, char** av) {
	ino_dir = (ac>1) ? av[1] : "/run/shm/ino-dir";
	fd_set rfs;
	int ino = inotify_init1(IN_CLOEXEC); FP1("ino: %d", ino);
	if (ino<0) perror("ino_ini"), exit(1);
	mkdir(ino_dir,0700); inotify_add_watch(ino, ino_dir, IN_CLOSE_WRITE|IN_CREATE|IN_DELETE);
	char buf[999];
	while(1) {
		FD_ZERO(&rfs); FD_SET(0, &rfs); FD_SET(ino, &rfs);
		int k = select(ino+1, &rfs, NULL, NULL, NULL); if (k<0) perror("sel"), exit(1);
		if (FD_ISSET(0, &rfs)) cmd(ino, buf, read(0, buf, 998));
		if (FD_ISSET(ino, &rfs)) hee(ino);
	}}
