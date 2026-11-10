
!$omp declare reduction(MaxPair:IndexPair:omp_out.Max(omp_in)) initializer(omp_priv=IndexPair(-1.e300,-1))
!$omp declare reduction(MinPair:IndexPair:omp_out.Min(omp_in)) initializer(omp_priv=IndexPair(1.e300,-1))
!$omp parallel for private(k,gx,ys) reduction(MaxPair:the_max) reduction(MinPair:the_min) schedule(static)
      do i = 1, this%nvecs
         ys = this%y(i)*this%stat(i)
         gx = this%grad(i)
         if (ys >= 0.0d0) then
            if (gx < gmin) then
               gmin = gx
               imin = i
            end if
         end if
         if ( ys <= 0.0d0) then
            if (gx > gmax) then
               gmax = gx
               imax = i
            end if
         end if
      end do
#else
