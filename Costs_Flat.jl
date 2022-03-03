

# x = Quadratmeter
# λ = Gewichtung des gemeinsamen Wohnraums

function cost(qm, gewichtung)
    Anteil = (qm + 90. /7 * gewichtung)/(sum(Zimmerqm) + 90. * gewichtung)
    round(Anteil * 2390., digits = 2)
end

Zimmerqm = [35.4, 24, 24, 18.8, 18, 15.65, 14.15]
costs0 = (cost.(Zimmerqm, 0))
costs1 = (cost.(Zimmerqm, 1))
costs2 = (cost.(Zimmerqm, 2))
costs15 = (cost.(Zimmerqm, 1.5))
costs25 = (cost.(Zimmerqm, 2.5))
costs3 = (cost.(Zimmerqm, 3))
costs4 = (cost.(Zimmerqm, 4))


total = [15, 40, 40, 45, 30, 30, 20] .+costs25

