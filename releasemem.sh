echo before:&&free -h && sync && echo 3 > /proc/sys/vm/drop_caches && echo after:&&free -h
