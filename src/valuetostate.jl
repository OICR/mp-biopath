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
        return "2"
    elseif downregulatedCutoff > value
        return "0"
    else
        return "1"
    end
end

end
