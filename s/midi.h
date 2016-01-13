#ifndef __qwe_midi_h__
#define __qwe_midi_h__

extern unsigned int midi_bv;
extern int midi_fd[32];
void midi_input(int i);
int midi_cmd(const char *s);

#endif // __qwe_midi_h__
