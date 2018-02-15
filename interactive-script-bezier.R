################################################################################################################################
#
# R script (for RStudio) for interactively visualize the Bezier curves corresponding to various discretized parameter values.
#
# Copyright (C) 2015  Dan Dediu
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
################################################################################################################################


################################################################################################################################
#
# This script *must* be run within RStudio (https://www.rstudio.com/)
#
################################################################################################################################


# Loading the 51-steps discretized parameter space of Bezier curves
# This can be *very* slow (about 20 mins on a Core i7 2630QM @2.0GHz 8Go RAM Ubuntu 14.04 64 bits with R 3.2.2 machine) and requires quite a bit of RAM (about 6Go on the same machine)!:
bezier.file <- xzfile("./generated-bezier.tsv.xz");
bezier <- read.table(bezier.file, header=TRUE, sep="\t", quote="");

# Extract the x-coordinates from the column names:
xs <- as.numeric(vapply(names(bezier)[6:ncol(bezier)], function(s){ substr(s,2,nchar(s)) }, character(1)));

# Optimization: cache the unique parameter values:
angle.values <- unique(bezier$angle); conc.values <- unique(bezier$conc); fronting.values <- unique(bezier$fronting); weigth.values <- unique(bezier$weigth); 

# Plot a given Bezier curve given either as a row in the grid.params data.frame or as the values of the four parameters:
.plot.Bezier.curve <- function( row=NA, # can be gives ad the row in the grid.params data.frame, or
                                angle=NA, conc=NA, fronting=NA, weigth=NA, # as the actual parameter values which uniquely select a row in the grid.params data.frame
                                tol=1e-3, # in which case we need a small tolerance when comparing floating point numbers for equality
                                rotate=NA, xy.ratio=NA # generic rotation and xy scaling for plotting
                              )
{
  if( is.na(row) && is.na(angle) && is.na(conc) && is.na(fronting) && is.na(weigth) ) return (FALSE);
  if( is.na(row) && (!is.na(angle) && !is.na(conc) && !is.na(fronting) && !is.na(weigth)) )
  {
    # Select the case using the closes actual parameter values to the ones given:
    angle.val    <- angle.values[    abs(angle.values -    angle)    <= tol ];
    conc.val     <- conc.values[     abs(conc.values -     conc)     <= tol ];
    fronting.val <- fronting.values[ abs(fronting.values - fronting) <= tol ];
    weigth.val   <- weigth.values[   abs(weigth.values -   weigth)   <= tol ];
    row <- which((bezier$angle == angle.val) & (bezier$conc == conc.val) & (bezier$fronting == fronting.val) & (bezier$weigth == weigth.val));
    if( length(row) != 1 )
    {
      # Nothing really to plot:
      plot( c(0,1), c(0,1), type="n", xlab="", ylab="", axes=FALSE,
            main=paste0("a=",sprintf("%.2f",angle)," c=",sprintf("%.2f",conc)," f=",sprintf("%.2f",fronting)," w=",sprintf("%.2f",weigth)));
      abline(v=seq(from=0, to=1, length.out=5), col=gray(0.8), lty="dashed");
      abline(h=seq(from=0, to=1, length.out=5), col=gray(0.8), lty="dashed");
      points( c(0,1), c(0,1), type="l", col="red"); points( c(0,1), c(1,0), type="l", col="red");
      mtext(c("back","front"), side=1, line=0, at=c(0,1), cex=0.75, col=gray(0.5));
      mtext(c("bottom","top"), side=2, line=0, at=c(0,1), cex=0.75, col=gray(0.5));
      return (FALSE);
    }
  }
  
  # The rotation and scaling:
  if( is.na(rotate) ) rotate <- -0.3217506; #  fixed rotation
  if( is.na(xy.ratio) ) xy.ratio <- 1/bezier$ratio[row];
  
  # The y coordinates:
  ys = as.numeric(bezier[row,-(1:5)]);
  
  # Rescale and rotate the plot:
  cos.rotate = cos(rotate); sin.rotate = sin(rotate); # cache for speedup
  M = matrix(c( cos.rotate, -sin.rotate,
                sin.rotate,  cos.rotate), byrow=TRUE, nrow=2);
  coords = matrix(c(xs-0.5, ys-0.5), byrow=TRUE, nrow=2);
  coords[2,] = coords[2,] * xy.ratio;
  coords2 = (M %*% coords);
  
  range.xs = range(coords2[1,]); range.ys = range(coords2[2,]); # speedups
  plot( coords2[1,], coords2[2,], type="n", xlab="", ylab="", axes=FALSE,
        main=paste0("a=",sprintf("%.2f",bezier$angle[row])," c=",sprintf("%.2f",bezier$conc[row])," f=",sprintf("%.2f",bezier$fronting[row])," w=",sprintf("%.2f",bezier$weigth[row])));
  abline(v=seq(from=range.xs[1], to=range.xs[2], length.out=5), col=gray(0.8), lty="dashed");
  abline(h=seq(from=range.ys[1], to=range.ys[2], length.out=5), col=gray(0.8), lty="dashed");
  points( coords2[1,], coords2[2,], type="l", col="blue");
  mtext(c("back","front"), side=1, line=0, at=range.xs, cex=0.75, col=gray(0.5));
  mtext(c("bottom","top"), side=2, line=0, at=range.ys, cex=0.75, col=gray(0.5));
        
  return (TRUE);
}
# TEST: .plot.Bezier.curve( row=123 )
# TEST: .plot.Bezier.curve( angle=0, conc=0, fronting=0.02, weigth=0.38 )

# Interactive exploration of the parameters:
library(manipulate);
manipulate( invisible(.plot.Bezier.curve(angle=angle, conc=conc, fronting=fronting, weigth=weigth)),
            angle=slider(0,1,initial=0.5,label="angle",step=0.02),
            conc=slider(0,1,initial=0.5,label="conc",step=0.02),
            fronting=slider(0,1,initial=0.5,label="fronting",step=0.02),
            weigth=slider(0,1,initial=0.5,label="weigth",step=0.02));
