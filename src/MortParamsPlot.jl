module MortParamsPlot

using DataArrays, DataFrames, DataFramesMeta, LaTeXStrings, LifeTable, PyCall, PyPlot, SQLite, LsqFit, GLM
import YAML

export LongTable, FitFrames, FitParams, FitParamsNonnorm, ParamsPlot, ObspredPlot, PredictMortalityPattern, PredictCauseLife, ObsCauseLife, ZeroAdd, CoeffForm

mgendir = expanduser("~/mortchartgen/")
conf = YAML.load_file(joinpath(mgendir, "chartgen.yaml"))
paramlabs = ["log_r0" =>  "\\mathrm{log}(r_0)", 
            "log_r0alpha" =>  "\\mathrm{log}\\frac{r0}{\\alpha}", 
            "alpha" =>  "\\alpha", "I(a - 1)" =>  "(a-1)", "a" =>  "a", 
            "trans_atau" =>  "\\mathrm{log}\\frac{a}{\\tau}-(a-1)\\mathrm{log}(\\tau)",
            "minalog_tau" =>  "-a\\mathrm{log}(\\tau)"]
PyDict(matplotlib["rcParams"])["axes.formatter.use_locale"] = true
matplotlib[:style][:use]("ggplot")

function PlotConf(func, functype, k)
	confdict = ["gompertz" =>  ["funclab" =>  "Gompertz", 
        "xcol" =>  ["rate" =>  "alpha", "surv" =>  "alpha"], 
        "ycol" =>  ["rate" =>  "log_r0", "surv" =>  "log_r0alpha"], "i" =>  -k],
        "weibull" =>  ["funclab" =>  "Weibull", 
	"xcol" =>  ["rate" =>  "I(a - 1)", "surv" =>  "a"],
        "ycol" =>  ["rate" =>  "trans_atau", "surv" =>  "minalog_tau"], "i" =>  exp(-k)],
        "functypelab" =>  ["rate" =>  "dödstal", "surv" =>  "överlevnad"]]
	return ["func" => ["funclab" => confdict[func]["funclab"], 
	"xcol" => confdict[func]["xcol"][functype],
	"ycol" => confdict[func]["ycol"][functype], "i" => confdict[func]["i"]],
	"functypelab" => confdict["functypelab"][functype]]
end

function LongTable(indbfile = joinpath(mgendir, "chartgen.db"))
	indb = SQLiteDB(indbfile)
	deaths_rs = query(indb, "select * from deaths")
	pop_rs = query(indb, "select * from pop")
	deaths_df = DataFrame(deaths_rs.values, convert(Vector{Symbol}, 
	deaths_rs.colnames))
	pop_df = DataFrame(pop_rs.values, convert(Vector{Symbol}, 
	pop_rs.colnames))
	agealias_arr = AgeAlias([2:25; "36sum"; "222sum"; "2325sum"])
	age_arr = [0:4; 5:5:95; 1; 0; 85]
	ages = DataFrame(Age = age_arr, Agealias = agealias_arr)
	deaths_long = stack(deaths_df, agealias_arr, [:Cause, :Sex, :Year, :Country])
	pop_long = stack(pop_df, agealias_arr, [:Sex, :Year, :Country])
	rename!(deaths_long, {:variable=>:Agealias, :value=>:Deaths})
	rename!(pop_long, {:variable=>:Agealias, :value=>:Pop})
	deaths_long = join(deaths_long, ages, on = :Agealias)
	pop_long = join(pop_long, ages, on = :Agealias)
	return [:deaths => deaths_long, :pop => pop_long]
end

function AgeAlias(suffices)
	convert(Vector{Symbol}, map((x)->"Pop$x", suffices))
end


function FitFrames(indict, sex, ages, ageformat, dimension, cause, country, year)

	alias = [0 => collect(2:25), 1 => [2:22; "2325sum"], 
	2 => [2; "36sum"; 7:22; "2325sum"]]
	alias_af = AgeAlias(alias[ageformat])

	function CommCond(frame)
		subframe = @ix(frame, (:Sex .== sex))
		subframe = @ix(subframe, findin(:Age, ages))
		subframe = @ix(subframe, findin(:Agealias, alias_af))
	end
	
	deaths = CommCond(indict[:deaths])
	pop = CommCond(indict[:pop])

	function DimSymb(col)
		convert(Vector{Symbol}, map((x)->"n$x", col))
	end

	if dimension == :Cause 
		deahts = @ix(deaths, findin(:Cause, cause))
		deaths = @ix(deaths , (:Country .== country) & (:Year .== year),
		[:Cause, :Deaths, :Age])
		pop = @ix(pop, (:Country .== country) & (:Year .== year), 
		[:Pop, :Age])
		deaths[:Cause] = DimSymb(deaths[:Cause])
		deaths = unstack(deaths, :Age, :Cause, :Deaths)
	elseif dimension == :Country
		deaths = @ix(deaths, findin(:Country, country))
		pop = @ix(pop, findin(:Country, country))
		deaths = @ix(deaths, (:Cause .== cause) & (:Year .== year),
		[:Country, :Deaths, :Age])
		pop =  @ix(pop, :Year .== year, [:Country, :Pop, :Age])
		deaths[:Country] = DimSymb(deaths[:Country])
		pop[:Country] = DimSymb(pop[:Country])
		deaths = unstack(deaths, :Age, :Country, :Deaths)
		pop = unstack(pop, :Age, :Country, :Pop)
	elseif dimension == :Year
		deaths = @ix(deaths, findin(:Year, year))
		pop = @ix(pop, findin(:Year, year))
		deaths = @ix(deaths, (:Country .== country) & (:Cause .== cause),
		[:Year, :Deaths, :Age])
		pop = @ix(pop, :Country .== country, [:Year, :Pop, :Age])
		deaths[:Year] = DimSymb(deaths[:Year])
		pop[:Year] = DimSymb(pop[:Year])
		deaths = unstack(deaths, :Age, :Year, :Deaths)
		pop = unstack(pop, :Age, :Year, :Pop)
	elseif dimension == :Age
		deaths = @ix(deaths, findin(:Year, year))
		pop = @ix(pop, findin(:Year, year))
		deaths = @ix(deaths, (:Country .== country) & (:Cause .== cause),
		[:Year, :Deaths, :Age])
		pop = @ix(pop, :Country .== country, [:Year, :Pop, :Age])
		deaths[:Age] = DimSymb(deaths[:Age])
		pop[:Age] = DimSymb(pop[:Age])
		deaths = unstack(deaths, :Year, :Age, :Deaths)
		pop = unstack(pop, :Year, :Age, :Pop)
	end

	return [:sex => sex, :ages => ages, :dimension => dimension, 
	:cause => cause, :country => country, :year => year, 
	:deaths => deaths, :pop => pop]

end

function FitParams(indict::Array, func, functype)
	indict_ca = indict[1]
	indict_tot = indict[2]
	cadeaths = indict_ca[:deaths]
	totdeaths = indict_tot[:deaths]
	pop = indict_tot[:pop]
	cols = size(cadeaths, 2)
	dimension = indict_ca[:dimension]
	sex = indict_ca[:sex]
	
	camsfs = Dict[]
	calts = DataFrame[]
	for i in 2:cols
		if(dimension == :Country || dimension == :Year)
			totcol = totdeaths[i]
			popcol = pop[i]
		elseif dimension == :Cause
			totcol = totdeaths[2]
			popcol = pop[2]
		end

		tot = DataFrame(Age = totdeaths[1], 
		Pop = popcol, Deaths = totcol)
		caprop = cadeaths[i]./totcol
		lt = PeriodLifeTable(tot, sex, false)
		calife = CauseLife(lt, caprop)
		calt = PeriodLifeTable(calife, sex, false, "rate")
		push!(calts, calt)
		camsf = MortSurvFit(calt, cadeaths[i], func, functype)
		push!(camsfs, camsf)
	end
	return [:indict => indict_ca, :lts => calts, :msfs => camsfs,
	:parafit => MortSurvParamsFit(camsfs), :func => func, :functype => functype]
end

function FitParams(indict::Dict, func, functype)
	deaths = indict[:deaths]
	pop = indict[:pop]
	cols = size(deaths, 2)
	dimension = indict[:dimension]
	sex = indict[:sex]
	
	msfs = Dict[]
	lts = DataFrame[]
	for i in 2:cols
		if(dimension == :Country || dimension == :Year)
			popcol = pop[i]
		elseif dimension == :Cause
			popcol = pop[2]
		end
		deathpop = DataFrame(Age = deaths[1], 
		Pop = popcol, Deaths = deaths[i])
		lt = PeriodLifeTable(deathpop, sex, false)
		push!(lts, lt)
		msf = MortSurvFit(lt, deaths[i], func, functype)
		push!(msfs, msf)
	end
	return [:indict => indict, :lts => lts, :msfs => msfs,
	:parafit => MortSurvParamsFit(msfs), :func => func, :functype => functype]
end

function PredictAgeRate(indict)
	deaths = indict[:deaths]
	pop = indict[:pop]
	cols = size(deaths, 2)
	model(x, p) = p[1] * exp(p[2].*x)

	pars = []
	obsrates = []
	for i in 2:cols
		year = convert(Array{Int}, deaths[1] - deaths[1][1]) 
		rate = convert(Array{Float64}, deaths[i]./pop[i])
		push!(obsrates, rate)
		fit = curve_fit(model, year, rate, [rate[1], log(rate[2])-log(rate[1])])
		push!(pars, fit.param)
	end
	return Dict(:model => model, :pars => pars, :obsrates => obsrates, :yroffset => deaths[1][1]) 
end

function PredictYear(preddict, predyear)
	preds = []
	for i in 1:size(preddict[:pars], 1)
		push!(preds, preddict[:model](predyear - preddict[:yroffset], preddict[:pars][i]))
	end
	return preds
end

function PredictTotal(indict, sex, ages, ageformat, country, obsyears, predyears)
	totdict =  FitFrames(indict, sex, ages, ageformat, :Age, "all", country, obsyears)
	cols = size(totdict[:deaths], 2)
	totpreddict = PredictAgeRate(totdict)
	rateframes = []
	for predyear in predyears
		rateframe = DataFrame(Age = ages, all = 0)
		totpred = PredictYear(totpreddict, predyear)
		rateframe[:all] = totpred
		push!(rateframes, rateframe)
	end
	return Dict(:sex => sex, :totpreddict => totpreddict, :rateframes => rateframes)
end

function PredictMortalityPattern(indict, sex, ages, ageformat, causes, country, obsyears, predyears, method = "curve_fit")
	othdict =  FitFrames(indict, sex, ages, ageformat, :Age, "all", country, obsyears)
	cols = size(othdict[:deaths], 2)
	if method == "curve_fit"
		predratefunc = PredictAgeRate
		predyrfunc = PredictYear
	elseif method == "svd"
		predratefunc = SvdCauseRate 
		predyrfunc = SvdPredictYear
	end
	capreddicts = Dict()
	for cause in causes
		cadict = FitFrames(indict, sex, ages, ageformat, :Age, cause, country, obsyears)
		capreddicts[cause] = predratefunc(cadict)
		for i in 2:cols
			othdict[:deaths][i] = othdict[:deaths][i] .- cadict[:deaths][i]
		end
	end
	capreddicts["oth"] = predratefunc(othdict)
	rateframes = []
	for predyear in predyears
		rateframe = DataFrame(Age = ages, all = 0)
		for cause in causes
			capred = predyrfunc(capreddicts[cause], predyear)
			rateframe[:all] = rateframe[:all] .+ capred
			casym = convert(Symbol, cause)
			rateframe[casym] = capred
		end
		othpred = predyrfunc(capreddicts["oth"], predyear)
		rateframe[:all] = rateframe[:all] .+ othpred
		rateframe[:oth] = othpred
		push!(rateframes, rateframe)
	end

	return Dict(:sex => sex, :capreddicts => capreddicts, :rateframes => rateframes)
end

function PredictCauseLife(mortpdict, causes)
	caframes = []
	ltframes = []
	for frame in mortpdict[:rateframes]
		carates = 0
		lt = PeriodLifeTable(frame, mortpdict[:sex], true, "rate")
		for cause in causes
			carates = carates .+ frame[convert(Symbol, cause)]
		end
		ca = CauseLife(lt, carates./frame[:all])
		push!(caframes, ca)
		push!(ltframes, lt)
	end
	return Dict(:caframes => caframes, :ltframes => ltframes)
end

function ObsCauseLife(indict, sex, ages, ageformat, cause, country, years)
	totdict = FitFrames(indict, sex, ages, ageformat, :Year, "all", country, years)
	cadict = FitFrames(indict, sex, ages, ageformat, :Year, cause, country, years)
	cols = size(totdict[:deaths], 2)
	caframes = []
	for i in 2:cols
		totframe = DataFrame(Age = totdict[:deaths][:Age], Pop = totdict[:pop][i],
		Deaths = totdict[:deaths][i])
		lt = PeriodLifeTable(totframe, sex)
		caframe = CauseLife(lt, cadict[:deaths][i] ./ totdict[:deaths][i])
		push!(caframes, caframe)
	end
	return Dict(:totdict => totdict, :cadict => cadict, :caframes => caframes)
end

function ZeroAdd(x, c = 0.5)
	if x == 0
		return c
	else
		return x
	end
end

function SvdCauseRate(indict)
	cols = size(indict[:deaths], 2)
	yroffset = indict[:year][1]
	deathsarr = map(ZeroAdd, transpose(convert(Array{Float64}, indict[:deaths][2:cols])))
	poparr = transpose(convert(Array{Float64}, indict[:pop][2:cols]))
	logratearr = log(deathsarr ./ poparr)
	rowmeans = []
	for i in 1:cols-1
		rowmean = mean(logratearr[i, :])
		push!(rowmeans, rowmean)
	end
	diffarr = convert(Array{Float64}, logratearr .- rowmeans)
	svddiff = svd(diffarr)
	svsq = svddiff[2].^2
	svfit = svsq ./ sum(svsq)
	ktframe = DataFrame(Years = indict[:year] .- yroffset, Kts = svddiff[3][:, 1])
	ktmodel = glm(Kts~Years, ktframe, Normal(), IdentityLink())
	return Dict(:rowmeans => rowmeans, :logratearr => logratearr, :svddiff => svddiff, :svfit => svfit, :ktmodel => ktmodel, :yroffset => yroffset)
end

function SvdPredictYear(indict, predyear, method = "glm")
	svddiff = indict[:svddiff]
	lastobsyear = size(svddiff[3], 1) + indict[:yroffset]
	if method == "glm"
		ktpred = coef(indict[:ktmodel])[1] + 
		coef(indict[:ktmodel])[2] * (predyear - indict[:yroffset])
	elseif method == "arima"
		ktpred = svddiff[3][end, 1] + (predyear - lastobsyear) * 
		((svddiff[3][end, 1] - svddiff[3][1, 1])/size(svddiff[3], 1))
	end
	bxpred = svddiff[1][:,1]
	logpreds = convert(Array{Float64}, indict[:rowmeans] .+ bxpred .* svddiff[2][1] .* ktpred)
	return exp(logpreds)
end

CoeffForm(coeff) = replace("$(round(coeff, 3))", ".", "{,}")

function RSquared(fit)
    sstot = sum((fit.model.rr.y - mean(fit.model.rr.y)).^2)
    ssres = sum((fit.model.rr.y - fit.model.rr.mu).^2)
    return (1-(ssres/sstot))
end

function ParamsPlot(fitdict, plot = "params")
	yrlab = "År"
	ages = fitdict[:indict][:ages]
	year = fitdict[:indict][:year]
	cause = fitdict[:indict][:cause]
	country = fitdict[:indict][:country]
	sex = fitdict[:indict][:sex]
	dimension = fitdict[:indict][:sex]
	causealias = conf["causes"][cause]["alias"] 
	sexalias = conf["sexes"][sex]["alias"] 
	countryalias = conf["countries"][country]["alias"] 
	trans_params = fitdict[:parafit]["trans_params"]
	font = ["size" => 14]
	coef_vec = coef(fitdict[:parafit]["fit"])
	b = coef_vec[1]
	k = coef_vec[2]
	rsquared = RSquared(fitdict[:parafit]["fit"])
	plotconf = PlotConf(fitdict[:func], fitdict[:functype], k)
	i = plotconf["func"]["i"]
	xcol = plotconf["func"]["xcol"]
	ycol = plotconf["func"]["ycol"]
	funclab = plotconf["func"]["funclab"]
	functypelab = plotconf["functypelab"]
	if dimension == :Cause
		dimstring = "per orsak"
	elseif dimension == :Country
		dimstring = "per befolkning"
	elseif dimension == :Year
		dimstring = "$(year[1])\u2013$(year[size(year, 1)])"
	end
	plottitle = string("$(funclab)analys $(functypelab) $(causealias)\n", 
	"$(sexalias) [$(ages[1]), $(ages[size(ages, 1)])] $(countryalias) ", 
	"$(dimstring)")
	paramstring = latexstring(paramlabs[ycol], "\\approx", CoeffForm(k),
	paramlabs[xcol], "+", CoeffForm(b))
	rsqstring = latexstring("R^2\\approx", CoeffForm(rsquared))
	istring = latexstring("i\\approx", CoeffForm(i))
	close()
	fig = figure()
	ax = fig[:add_subplot](111)
	if plot == "params"
		ax[:scatter](trans_params[:X], trans_params[:Y])
		ax[:set_xlabel](latexstring(paramlabs[xcol]))
		ax[:set_ylabel](latexstring(paramlabs[ycol]))
		for (string, ypos) in [(paramstring, 0.9), (rsqstring, 0.84), (istring, 0.78)]
			ax[:text](0.98, ypos, string, transform=ax[:transAxes], 
			fontdict = font, ha = "right")
		end
	elseif plot == "trend"
		ax[:plot](year, trans_params[:X])
		ax[:set_xlabel](yrlab)
		ax[:set_ylabel](latexstring(paramlabs[xcol]))
		ax[:set_xlim](year[1], year[(size(year, 1))])
	end
	ax[:set_title](plottitle)
end

function ObspredPlot(fitdict, points, ages, 
	trans = "none", trans_time_coords = false, trans_y_coords = false)
	functype = fitdict[:functype]
	year = fitdict[:indict][:year]
	cause = fitdict[:indict][:cause]
	country = fitdict[:indict][:country]
	sex = fitdict[:indict][:sex]
	causealias = conf["causes"][cause]["alias"] 
	sexalias = conf["sexes"][sex]["alias"] 
	countryalias = conf["countries"][country]["alias"] 
	plotconf = PlotConf(fitdict[:func], functype, 0)
	funclab = plotconf["func"]["funclab"]
	functypelab = plotconf["functypelab"]
	xlinlab = L"\mathrm{log}(t)"
	ratelinlab = L"\mathrm{log}(r(t))"
	survlinlab = L"\mathrm{log}[-\mathrm{log}[S(t)]]"
	plotlabs = ["none" => ["x" => L"t", "rate" => L"r(t)", "surv" => L"S(t)"],
        "gomp_lin" => ["x" => L"t", "rate" => ratelinlab, "surv" => survlinlab],
        "weib_lin" => ["x" => xlinlab,"rate" => ratelinlab, "surv" => survlinlab]]
	close()

	if functype == "rate"
		obscol = :m
	elseif functype == "surv"
		obscol = :l
	end

	fig = figure()
	ax = fig[:add_subplot](111)
	preddicts = Dict[]
	for point in points
		pointsym = symbol("n$point")
		pointind = findin(names(fitdict[:indict][:deaths]), [pointsym])[1] - 1
		obs_t = transdict(fitdict[:indict][:ages], 
		fitdict[:lts][pointind][obscol])
		pred_t = transdict(ages, Predict(ages, fitdict[:msfs][pointind]))
		push!(preddicts, pred_t)
		obsplot = ax[:plot](obs_t[trans]["x"], 
		obs_t[trans][functype], "o", label = string(point))
		curcolor = obsplot[1][:get_color]()
		predplot = ax[:plot](pred_t[trans]["x"], 
		pred_t[trans][functype], color = curcolor)
	end
	plottitle = string("Observerad vs förutsedd $(funclab)analys $(functypelab)\n", 
	"$(causealias) $(sexalias) $(countryalias)")
	legend(loc = 2, framealpha = 0.5)
	ax[:set_title](plottitle)

	if (trans_time_coords)
		ax[:xaxis][:set_major_locator](plt[:FixedLocator](preddicts[1][trans]["x"]))
		locs, labels = xticks()
		xticks(locs, preddicts[1]["none"]["x"])
		ax[:set_xlabel](plotlabs["none"]["x"])
		xticks(rotation = 90)
	else
		ax[:set_xlabel](plotlabs[trans]["x"])
	end
	if (trans_y_coords)
		ax[:yaxis][:set_major_locator](plt[:FixedLocator](preddicts[1][trans][functype]))
		locs, labels = yticks()
		trans_y_ticks = map((x)->replace(string(round(x, 6)), ".", ","), 
		preddicts[1]["none"][functype])
		yticks(locs, trans_y_ticks)
		ax[:set_ylabel](plotlabs["none"][functype])
	else
		ax[:set_ylabel](plotlabs[trans][functype])
	end
end

function transdict(x, y)
	xlin = log(x)
	ratelin = log(y)
	survlin = log(-log(y))
	return ["none" => ["x" => x, "rate" => y, "surv" => y],
		"gomp_lin" => ["x" => x, "rate" => ratelin, "surv" => survlin],
		"weib_lin" => ["x" => xlin, "rate" => ratelin, "surv" => survlin]]
end

end  #module
