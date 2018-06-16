

    # KRP W

    KRP_SubSampling_W = matrix(0,4,11)
	KRP_SubSampling_W[1,] = SubSamplingMatrixTest( KRP_5_W_rpi, 17*19 )
	KRP_SubSampling_W[2,] = SubSamplingMatrixTest( KRP_5_W_rpi, 17*19 )
	KRP_SubSampling_W[3,] = SubSamplingMatrixTest( KRP_5_W_rpi, 17*19 )
	KRP_SubSampling_W[4,] = SubSamplingMatrixTest( KRP_6_W25_rpi, 17*19 )

    # KRP WH
    KRP_SubSampling_WH = matrix(0,4,11)
	KRP_SubSampling_WH[1,] = SubSamplingMatrixTest( KRP_5_WH4AT1_rpi, 17*19 )
	KRP_SubSampling_WH[2,] = SubSamplingMatrixTest( KRP_5_WH4AT2_rpi, 17*19 )
	KRP_SubSampling_WH[3,] = SubSamplingMatrixTest( KRP_5_WH2AT_rpi, 17*19 )
	KRP_SubSampling_WH[4,] = SubSamplingMatrixTest( KRP_6_WH25_rpi, 17*19 )

    X71_SubSampling_WH = matrix(0,9,11)
    X71_SubSampling_WH[1,] = SubSamplingMatrixTest( bm2_A ) # * work
    X71_SubSampling_WH[2,] = SubSamplingMatrixTest( bm2_Q ) # * work
    X71_SubSampling_WH[3,] = SubSamplingMatrixTest( bm3_A ) # * work
    X71_SubSampling_WH[4,] = SubSamplingMatrixTest( bm3_Q ) # * work
    X71_SubSampling_WH[5,] = SubSamplingMatrixTest( bm3_QL ) # * work
    X71_SubSampling_WH[6,] = SubSamplingMatrixTest( bm4_S4A ) # * work
    X71_SubSampling_WH[7,] = SubSamplingMatrixTest( bm4_SQ ) # * work
    X71_SubSampling_WH[8,] = SubSamplingMatrixTest( bm6_SQ ) # * work
    X71_SubSampling_WH[9,] = SubSamplingMatrixTest( bm7_SQ ) # * work

    
    X71_SubSampling_W = matrix(0,9,11)
    X71_SubSampling_W[1,] = SubSamplingMatrixTest( bm2_W ) # * work
    X71_SubSampling_W[2,] = SubSamplingMatrixTest( bm2_W ) # * work
    X71_SubSampling_W[3,] = SubSamplingMatrixTest( bm3_W ) # * work
    X71_SubSampling_W[4,] = SubSamplingMatrixTest( bm3_W ) # * work
    X71_SubSampling_W[5,] = SubSamplingMatrixTest( bm3_W ) # * work
    X71_SubSampling_W[6,] = SubSamplingMatrixTest( bm4_SW ) # * work
    X71_SubSampling_W[7,] = SubSamplingMatrixTest( bm4_SW ) # * work
    X71_SubSampling_W[8,] = SubSamplingMatrixTest( bm6_SW ) # * work
    X71_SubSampling_W[9,] = SubSamplingMatrixTest( bm7_SW ) # * work


    X163_SubSampling_WH = matrix(0,7,11)
    X163_SubSampling_WH[1,] = SubSamplingMatrixTest( bP170_2_Q, 173*173 )
    X163_SubSampling_WH[2,] = SubSamplingMatrixTest( bP170_3_Q, 173*173 )
    X163_SubSampling_WH[3,] = SubSamplingMatrixTest( bP170_4_S4A, 173*173 )
    X163_SubSampling_WH[4,] = SubSamplingMatrixTest( bP170_4_SQ, 173*173 )
    X163_SubSampling_WH[5,] = SubSamplingMatrixTest( bP170_5_SQ, 173*173 )
    X163_SubSampling_WH[6,] = SubSamplingMatrixTest( bP170_6_SA1, 173*173 )
    X163_SubSampling_WH[7,] = SubSamplingMatrixTest( bP170_6_SA2, 173*173 )


    X163_SubSampling_W = matrix(0,7,11)
    X163_SubSampling_W[1,] = SubSamplingMatrixTest( bP170_2_W, 173*173 )
    X163_SubSampling_W[2,] = SubSamplingMatrixTest( bP170_3_W, 173*173 )
    X163_SubSampling_W[3,] = SubSamplingMatrixTest( bP170_4_SW, 173*173 )
    X163_SubSampling_W[4,] = SubSamplingMatrixTest( bP170_4_SW, 173*173 )
    X163_SubSampling_W[5,] = SubSamplingMatrixTest( bP170_5_SW, 173*173 )
    X163_SubSampling_W[6,] = SubSamplingMatrixTest( bP170_6_SW1, 173*173 )
    X163_SubSampling_W[7,] = SubSamplingMatrixTest( bP170_6_SW2, 173*173 )




	par(mfrow=c(2,3))

    # X71 W
    x = c(0.1,0.5,1,2,4,8,16,32,64) #,128,256)
	y = (colSums(X71_SubSampling_W)/9)[1:9]
	std = (colSds(X71_SubSampling_W))[1:9]

	plot(x,y,ylim=c(0.3,1.0))
	segments(x, y-std,x, y+std, col = "dark blue")
	epsilon = 1
	segments(x-epsilon,y-std,x+epsilon,y-std, col = "dark blue")
	segments(x-epsilon,y+std,x+epsilon,y+std, col = "dark blue")


    # X163 W
    x = c(0.1,0.5,1,2,4,8,16,32,64) #,128,256)
	y = (colSums(X163_SubSampling_W)/7)[1:9]
	std = (colSds(X163_SubSampling_W))[1:9]

	plot(x,y,ylim=c(0.3,1.0))
	segments(x, y-std,x, y+std, col = "dark blue")
	epsilon = 1
	segments(x-epsilon,y-std,x+epsilon,y-std, col = "dark blue")
	segments(x-epsilon,y+std,x+epsilon,y+std, col = "dark blue")

    # # KRP W
	x = c(0.1,0.5,1,2,4,8,16,32,64) #,128,256)
	y = (colSums(KRP_SubSampling_W)/4)[1:9]
	std = (colSds(KRP_SubSampling_W))[1:9]

	plot(x,y,ylim=c(0.3,1.0))
	segments(x, y-std,x, y+std, col = "dark blue")
	epsilon = 1
	segments(x-epsilon,y-std,x+epsilon,y-std, col = "dark blue")
	segments(x-epsilon,y+std,x+epsilon,y+std, col = "dark blue")


    # X71 WH
    x = c(0.1,0.5,1,2,4,8,16,32,64) #,128,256)
	y = (colSums(X71_SubSampling_WH)/9)[1:9]
	std = (colSds(X71_SubSampling_WH))[1:9]

	plot(x,y,ylim=c(0.6,1.0))
	segments(x, y-std,x, y+std, col = "dark blue")
	epsilon = 1
	segments(x-epsilon,y-std,x+epsilon,y-std, col = "dark blue")
	segments(x-epsilon,y+std,x+epsilon,y+std, col = "dark blue")


    # X163 WH
    x = c(0.1,0.5,1,2,4,8,16,32,64) #,128,256)
	y = (colSums(X163_SubSampling_WH)/7)[1:9]
	std = (colSds(X163_SubSampling_WH))[1:9]

	plot(x,y,ylim=c(0.6,1.0))
	segments(x, y-std,x, y+std, col = "dark blue")
	epsilon = 1
	segments(x-epsilon,y-std,x+epsilon,y-std, col = "dark blue")
	segments(x-epsilon,y+std,x+epsilon,y+std, col = "dark blue")

    # # KRP W

	x = c(0.1,0.5,1,2,4,8,16,32,64) #,128,256)
	y = (colSums(KRP_SubSampling_WH)/4)[1:9]
	std = (colSds(KRP_SubSampling_WH))[1:9]

	plot(x,y,ylim=c(0.6,1.0))
	segments(x, y-std,x, y+std, col = "dark blue")
	epsilon = 1
	segments(x-epsilon,y-std,x+epsilon,y-std, col = "dark blue")
	segments(x-epsilon,y+std,x+epsilon,y+std, col = "dark blue")

    