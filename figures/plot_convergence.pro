pro plot_convergence,tims,min_values,title=title

; routine to plot the convergence of the penalty function to
;   the optimum value


list=where(tims ne 0)
plot,tims[list]/60.,min_values[list],xtitle= $
     'Time (min)',ytitle= $
     'Penalty Function at Best Vertex (arb.)',title=title, $
     background=255,color=0,charsize=1.5
     
     
end