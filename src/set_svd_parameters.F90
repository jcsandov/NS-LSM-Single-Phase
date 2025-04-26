subroutine set_svd_parameters( mms , nms )

  integer, intent(in) :: mms , nms

  lda = max( 1 , mms )
  ldb = max( 1 , max( mms , nms ) )

  nlvl = max(0,int( log( real( min(mms,nms) , kind=rdf )/real(smlsiz+1,kind=rdf))/log(two) )+1)

  if( mms >= nms ) then      
    lwork  = 12*nms+2*nms*smlsiz+8*nms*nlvl+ nms*nrhs+(smlsiz+1)*(smlsiz+1)
  else
    lwork  = 12*mms+2*mms*smlsiz+8*mms*nlvl+ mms*nrhs+(smlsiz+1)*(smlsiz+1)
  end if
  
  liwork = max(1,3*min(mms,nms)*nlvl+11*min(mms,nms))
  
  allocate( work(  max( 1 , lwork  ) ) )
  allocate( iwork( max( 1 , liwork ) ) )
  allocate( s( max( 1 , min( mms , nms ) ) ) )

end subroutine set_svd_parameters