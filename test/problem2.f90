program problem2
implicit none
  double precision a,b,h,trap,mid,simp !a,bは始点と終点．今回は0と1にする．
  double precision,pointer::T(:,:)  !subroutine"Romberg"からの指示先引継ぎ用のポインタ
  integer,parameter::Jmax=3   !Romberg積分のJの最大値
  double precision,target:: Tlist(0:Jmax,0:Jmax)  !実際の積分値が入る配列
  integer N,k,i,j
  double precision,parameter::pi=4.0d0*(4.0d0*atan(1.0d0/5.0d0)-atan(1.0d0/239d0))
                                  !πの定義(Machin's formula)

  a=0.0d0
  b=1.0d0
  open(10,file='trap.dat',status='replace')
  open(11,file='mid.dat',status='replace')
  open(12,file='simp.dat',status='replace')
  do k=0,7
    N=2**k  !分割数
    h=1.0d0/N !刻み幅
    call trapezoid(a,b,h,N,trap)  !台形測
    call midpoint(a,b,h,N,mid) !中点則
    call simpson(a,b,h,N,simp)  ! Simpson則
    write(10,*) k,N,trap,abs(trap-pi/4.0d0)
    write(11,*) k,N,mid,abs(mid-pi/4.0d0)
    write(12,*) k,N,simp,abs(simp-pi/4.0d0)
  end do
  close(10)
  close(11)
  close(12)
  
  T=>Tlist      !TlistにTを代入
  call Romberg(Jmax,T)    !引数のTはポインタの指示先を引き渡す
  open(14,file='romberg.dat',status='replace')
  do i=0,Jmax   
    do j=0,i
      write(14,*) Tlist(i,j),abs(Tlist(i,j)-(pi/4.0d0))
    end do
    write(14,*)
  end do
  write(*,*) Tlist
contains

!-------------------------------------------------------------------------------------------

 double precision function func(x)  !被積分関数
 implicit none
  double precision x
    func=1.0d0/(1.0d0+x**2)
 end function

!-------------------------------------------------------------------------------------------

 subroutine midpoint(a,b,h,N,mid)  !複合中点則による計算
 implicit none
  double precision,intent(in)::a,b,h  !入力
  integer,intent(in)::N   !入力
  double precision,intent(out)::mid  !出力
  double precision sum    !local変数
  integer i
  
  sum=0.0d0
  do i=1,N
    sum=sum+func(a+h*(i-0.5d0))
  end do
  mid=h*sum
 end subroutine

!-----------------------------------------------------------------------------------------------

 subroutine trapezoid(a,b,h,N,trap) !複合台形則による計算
 implicit none
  double precision,intent(in)::a,b,h  !入力
  integer,intent(in)::N   !入力
  double precision,intent(out)::trap  !出力
  double precision sum    !local変数
  integer i
  
  sum=0.0d0
  do i=1,N-1
    sum=sum+func(a+i*h)
  end do
  trap=h*((func(a)+func(b))/2.0d0+sum)
 end subroutine

 subroutine simpson(a,b,h,N,simp) !複合Simpson則による計算
 implicit none
  double precision,intent(in)::a,b,h  !入力
  integer,intent(in)::N   !入力
  double precision,intent(out)::simp  !出力
  double precision sum1,sum2    !local変数
  integer i
  
  sum1=0.0d0
  sum2=0.0d0
  do i=1,N
    sum1=sum1+func(a+h*(i-0.5d0))
  end do
  do i=1,N-1
    sum2=sum2+func(a+i*h)
  end do
  
  simp=h/6.0d0*(func(a)+func(b)+4.0d0*sum1+2.0d0*sum2)
 end subroutine

!--------------------------------------------------------------------------------------------------

 double precision function Richard(T1,T2,j) !リチャードソン補外
 implicit none
    double precision T1,T2
    integer j
    
    Richard=(T1-4.0d0**(j)*T2)/(1.0d0-4.0d0**(j))
 end function

!--------------------------------------------------------------------------------------------------

 subroutine Romberg(Jmax,S)     !Romberg積分用の表をつくるsubroutine．表を示すポインタを引数とする．
 implicit none
  double precision,pointer::S(:,:)  !表を保存するための正方行列を示すポインタ
  integer,intent(in)::Jmax
  integer i,j,N
  double precision h,M
        !SはJmax+1行Jmax+1列の正方行列．ただし下三角分しか使わない．
  
  do i=0,Jmax   !表の一列目を作る．中点則と台形則より分割が倍の台形則を計算．
    N=2**i
    h=1.0d0/N
    if (i==0) then
      call trapezoid(0.0d0,1.0d0,h,N,S(i,i))
    else
      call midpoint(0.0d0,1.0d0,2.0d0*h,N/2,M)
      S(i,0)=(S(i-1,0)+M)/2.0d0
    end if
  end do
  
  do i=1,Jmax   !表を階段順に埋める
    do j=1,i
      S(i,j)=Richard(S(i-1,j-1),S(i,j-1),j)
    end do
  end do
 end subroutine
!-----------------------------------------------------------------------------------------------------
 end program

