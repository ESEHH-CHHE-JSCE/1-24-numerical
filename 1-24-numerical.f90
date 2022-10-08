!++++++++++++++++++++++++++++++++++++++++
! 1次元移流分散方程式の数値解析（有限差分法） 
!++++++++++++++++++++++++++++++++++++++++
program reidai082
      implicit none
      double precision :: c(1000),co(1000)
      double precision :: xmax,dx,dx2,xo,vel,al,tau,dm,dd,bc,tend,dt,tjikan,time1
      double precision :: o,p,q,r,s
      integer :: imin,imax,ii1,nmax,iout1,iout,i,n

! 出力ファイルの指定
      open(11,file='output.csv')

! 解析に用いる諸数値
      xmax = 30.0d0                              ! 解析領域の長さ (m)
      dx = 1.0d0                                 ! 差分格子間隔 (m)
      dx2 = dx ** 2                              ! 差分格子間隔の二乗
      imin = 1                                   ! 上流境界の格子点番号
      imax = xmax / dx + 1                       ! 下流境界の格子点番号

      xo = 10.0d0                                ! 濃度観測地点（m）
      ii1 = xo / dx + 1                          ! 濃度観測地点の格子点番号

      vel = 0.1d0                                ! 実流速 (m/day)
      al = 1.0d0                                 ! 縦分散長 (m)
      tau = 1.0d0                                ! 屈曲率（-） 
      dm = 1.0d-5                                ! 分子拡散係数（cm2/s）
      dd = al * abs(vel) + tau * dm *8.64d0      ! 分散係数 (m2/day)

      bc = 1.0d0                                 ! 上流境界の濃度 (mg/L)

      tend = 300.0                               ! 計算時間（day）
      dt = 1.0d0                                 ! 差分時間間隔（時間ステップ）（day）
      nmax = tend / dt                           ! 計算時間のステップ数

      tjikan = 1.0d0                             ! 結果をファイルに出力する時間間隔（day）
      iout1 = tjikan / dt                        ! 結果をファイルに出力する時間間隔のステップ数

! 初期条件
      time1 = 0.0d0
      iout = 0
      do i = imin, imax
       c(i) = 0.0d0
       if (i == 1) c(i) = bc
       co(i) = c(i)
      end do

! 濃度観測地点における初期時間と初期濃度のファイル出力
      write(11,*) 'elapsed-time(days), concentration(mg/L)'
      write(11,119) time1, c(ii1)
 119  format(f7.1,',',f11.7)

!///// 計算開始 /////
      do n = 1, nmax                            ! 時間ステップの繰り返し
       time1 = time1 + dt
       iout = iout + 1
        if (iout == iout1) then           
         write(*,120) time1
 120     format(f6.1,2x,'days')
        else
        end if

! 境界条件
       co(imin) = bc                            ! 上流境界条件
       co(imax) = 0.0d0                         ! 下流境界条件

! 上下流境界以外の濃度の計算
       do i = imin+1, imax-1
        o = co(i+1)
        p = co(i)
        q = co(i-1)
        r = dd * (o - 2 * p + q) * dt / dx2
        s = vel * (p - q) * dt / dx
        c(i) = p + r - s
       end do

! 既知濃度の更新
       do i = imin+1, imax-1
        co(i) = c(i)
       end do

! 一定時間ごとの濃度観測地点における時間と濃度のファイル出力
       if (iout == iout1) then
        write(11,119) time1, c(ii1)
        iout = 0
       else
       end if

      end do
!///// 計算終了 /////*

! 出力ファイルを閉じる
      close(11)
      close(12)

      stop
end program reidai082
