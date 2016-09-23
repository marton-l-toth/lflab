define trkdump
set $p = $arg0->g0
set $pp = 0
set $z = $arg0->g1
printf "%p: g0=%p g1=%p\n", $arg0, $p, $z
while $p && $pp != $z
	set $q = $p->m_u24.t
	set $c = $q.ty ? $q.ty : 48
	printf "%05x %c %p< %p >%p i%d j%d\n",$p->m_id,$c,$q.pv,$p,$p->m_next,$q.i,$q.j
	set $pp = $p
	set $p = $p->m_next
end
end
