LF_TMPROOT="$(readlink -f /dev/shm)"
[[ -d "$LF_TMPROOT" ]] || LF_TMPROOT="$(readlink -f /run/shm)"
[[ -d "$LF_TMPROOT" ]] || LF_TMPROOT="$(readlink -f /tmp)"
