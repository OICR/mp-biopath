module ValueToState

function getState(value_str, downregulatedCutoff, upregulatedCutoff)
    value = typeof(value_str) == Float64? value_str: parse(Float64, value_str)
    if upregulatedCutoff < value
        return "Up Regulated"
    elseif downregulatedCutoff > value
        return "Down Regulated"
    else
        return "Normal"
    end
end

function getStateNumber(value_str, downregulatedCutoff, upregulatedCutoff)
    value = parse(Float64, value_str)
    if upregulatedCutoff < value
        return "3"
    elseif downregulatedCutoff > value
        return "1"
    else
        return "2"
    end
end

end
