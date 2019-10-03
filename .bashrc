#PS1='[\[\e[33m\]#\##\[\033[01;04;32m\]\u@\[\033[32m\]\H:\[\e[36m\]\w\[\e[31m\]\$ \[\033[00m\]'
PS1='\[\e[33m\]\u\[\e[32m\]@\[\e[34m\]\h:\[\e[36m\]\W\[\e[01;031m\]\$\[\e[00m\]'
alias ll='ls -alF'
alias hisgrep='history | grep'
alias lstail='ls -ltrh | tail'
alias lshead='ls -ltrh | head'
alias la='ls -a'
alias rm='rm -i'
alias t20='tail -n 20'
alias ls="ls --color=auto"
alias l='ls -CF'
alias ..="cd .."
export HISTTIMEFORMAT='%F %T '
export HISTSIZE="5000"
export HISTFILESIZE=3000
export HISTFILE=~/$USER-$UID.log
shopt -s histappend
export PROMPT_COMMAND="history -a"
extract() {
    if [ -f $1 ] ; then
      case $1 in
        *.tar.bz2)   tar xjf $1     ;;
        *.tar.gz)    tar xzf $1     ;;
        *.bz2)       bunzip2 $1     ;;
        *.rar)       unrar e $1     ;;
        *.gz)        gunzip $1      ;;
        *.tar)       tar xf $1      ;;
        *.tbz2)      tar xjf $1     ;;
        *.tgz)       tar xzf $1     ;;
        *.zip)       unzip $1       ;;
        *.Z)         uncompress $1  ;;
        *.7z)        7z x $1        ;;
        *)     echo "'$1' cannot be extracted via extract()" ;;
         esac
     else
         echo "'$1' is not a valid file"
     fi
}
cl() { cd "$1"; ls;}
#PASS=$(cat ~/PASSWORD)
#alias su='echo ${PASS} > ~/pass; sudo -i -S < ~/pass 2> /dev/null; sudo -i -S'
#alias peregrine="~/SSHresponse peregrine.hpc.rug.nl p277903 Xmm205098"
#alias md-fw="~/sshHostWithUsernameAndPassword md-fw.housing.rug.nl minmin xmm205098"
#alias node10='~/SSHresponse 192.168.1.10 xuemm 111111'
#alias node20='~/SSHresponse 192.168.1.20 xuemm 111111'
#alias node30='~/SSHresponse 192.168.1.30 xuemm 111111'
#alias node40='~/SSHresponse 192.168.1.40 xuemm 111111'
#alias node88='~/SSHresponse 192.168.1.88 xuemm 111111'
#alias node99='~/SSHresponse 192.168.1.99 xuemm 111111'
#alias node88='~/SSHresponse 192.168.1.88 xuemm 111111'
#alias node99='~/SSHresponse 192.168.1.99 xuemm 111111'


#alias g2016='source /usr/local/gromacs-2016.1/bin/GMXRC'
alias g514='source /home/xuemm/software/gromacs-5.1.4/bin/GMXRC'
alias g507='source /home/xuemm/software/gromacs-5.0.7/bin/GMXRC'
alias g456='source /home/xuemm/software/gromacs-4.5.6/bin/GMXRC'
alias g467='source /home/xuemm/software/gromacs-4.6.7/bin/GMXRC'
#export MEMDIR=/coarse/minmin/membrane
export PATH=/opt/TurboVNC/bin/:/home/xuemm/bin:/public/software/compiler/intel/2017/bin/:$PATH
export LD_LIBRARY_PATH=/public/shared/:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/public/software/intel/composer_xe_2015.0.090/compiler/lib/intel64/:/public/software/intel/composer_xe_2015.0.090/mkl/lib/intel64:/usr/lib/:/usr/lib64/:/public/software/vmd/:/public/home/xuemm/software/gromacs-5.0.7/build/lib/:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/usr/lib/:/usr/lib64/:/public/software/vmd/:/home/xuemm/software/gromacs-5.0.7/build/lib/:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/usr/lib/:/usr/lib64/:/public/software/vmd/:/home/xuemm/software/gromacs-5.0.7/build/lib/:/public/software/mpi/openmpi1.6.5-intel/lib:/public/software/intel/composer_xe_2015.0.090/compiler/lib/intel64:/public/software/intel/composer_xe_2015.0.090/mpirt/lib/intel64:/public/software/intel/composer_xe_2015.0.090/ipp/../compiler/lib/intel64:/public/software/intel/composer_xe_2015.0.090/ipp/lib/intel64:/public/software/intel/composer_xe_2015.0.090/compiler/lib/intel64:/public/software/intel/composer_xe_2015.0.090/mkl/lib/intel64:/public/software/intel/composer_xe_2015.0.090/tbb/lib/intel64/gcc4.4:/usr/local/cuda-6.5/lib64:/usr/local/netpbm/lib:
export LESSOPEN="| /usr/bin/src-hilite-lesspipe.sh %s"
#export LESSOPEN="|~/.lessfilter %s"
export LESS=" -R "
alias cp='cp -i'
alias mv='mv -i'
alias vi='vim'
#alias cdmd88home='cd /md88-home/MD88/home'
#alias cdcoarse='cd /coarse/MD88/coarse'
#alias membranedirectory='cd /coarse/minmin/membrane/ && ls -d */'
#source /etc/rc.local
alias ovito='vglrun /home/xuemm/software/ovito-2.9.0-x86_64/bin/ovito'
