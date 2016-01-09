#ifndef __qwe_contrib_h__
#define __qwe_contrib_h__

class ANode;
typedef void (c_inifun)(ANode*);
struct contrib_ent { const char *s; c_inifun* f; };
extern contrib_ent contrib_list[];

#endif // __qwe_contrib_h__
