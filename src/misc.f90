module misc
  !$ use omp_lib
  implicit none
  public
  ! private
  ! public :: start_clock, stop_clock, assert

  real(8),private :: time ! [s]
  
contains
  subroutine start_clock
    call cpu_time(time)
    !$ time = omp_get_wtime()
  end subroutine start_clock

  subroutine stop_clock(t)
    real(8) :: t
    
    call cpu_time(t)
    !$ t = omp_get_wtime() 

    t = t - time
  end subroutine stop_clock

  subroutine assert(flag, message)
    logical,intent(in) :: flag
    character(*),intent(in),optional :: message

    if(.not. flag) then
       write(*,*) "Assertion failed"
       if(present(message)) write(*,*) message
         
#ifdef IFORT
       call tracebackqq
#endif
       
       call abort
    
    end if
  end subroutine assert

  !> 昇順クイックソート
  !!\param[in] n 配列サイズ
  !!\param[inout] x ソートする配列
  recursive subroutine quicksort(n, x)
    integer,intent(in) :: n
    integer,intent(inout) :: x(n)
    integer :: p,i,imin,imax,it,is
    if(n == 1) return    !配列長さ1ならそのまま返す、再帰の終点
    p = x(1)   !一個目を基準値に
    imin = 1   !前からスキャンしている場所
    imax = n   !後ろからスキャンしている場所
    do while(imin < imax)
       !スキャンと交換をするループ
       ! iminがimaxより前にある間だけで大丈夫、重なったら終わり
       do while(x(imin) < p .and. imin < imax)
          !前からみていって、p未満なら次の場所に進む、p以上ならとどまり交換を待つ
          ! 最初は必ずimin = 1で止まり、交換の対象となる
          imin = imin + 1
       enddo
       do while(p <= x(imax) .and. imin < imax)
          !後ろから見ていって、p以上なら次の場所に進む、p未満ならとどまり交換を待つ
          imax = imax - 1
       enddo
       !前からスキャンと後ろからスキャンが両方止まったら交換の合図、両者の値を入れ替える
       if(imin < imax)then
          it = x(imax)
          x(imax) = x(imin)
          x(imin) = it
       endif
    enddo
    ! 前からスキャンと後ろからスキャンの位置が重なっていたらスキャン終わり、配列を分割
    ! imax (この時点ではiminと同じはず)で分割するが、ちょうどimaxをどちらの配列に入れるかはお好みで、今回は後ろに
    is = imax - 1
    if(is == 0) is = 1
    ! 基準値より小さいほうをさらにソート、再帰呼び出し
    call quicksort(is, x(1:is))
    ! 基準値以上のほうをさらにソート、再帰呼び出し
    call quicksort(n - is, x(is + 1:n))
  end subroutine quicksort
  
end module misc
