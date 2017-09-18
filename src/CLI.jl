# Copyright (c) 2015 Michael Eastwood
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

__precompile__()

module CLI

export Command, Option

immutable Command
    name::UTF8String
    help::UTF8String
end

immutable Option
    flag::UTF8String
    help::UTF8String
    T::Type
    required::Bool
    min::Number # The minimum number of arguments passed to this option
    max::Number # The maximum number of arguments passed to this option
    conflicts::Vector{UTF8String} # A list of options that this one conflicts with
end

Option(flag,help;
       T=Void,
       required=false,
       min=0,
       max=0,
       conflicts=UTF8String[]) = Option(flag,help,T,required,min,max,conflicts)

function Base.print(io::IO,command::Command)
    print(io,"  ")
    print(io,rpad(command.name,16," "))
    print(io,replace(command.help,"\n","\n"*" "^18))
end

function Base.print(io::IO,option::Option)
    print(io,"  ")
    print(io,rpad(option.flag,16," "))
    print(io,replace(option.help,"\n","\n"*" "^18))
    option.required && print(io," (required)")
end

name = ""
banner = ""
commands = Command[]
options = Dict{AbstractString,Vector{Option}}()

set_name(str::AbstractString) = (global name; name = str)
set_banner(str::AbstractString) = (global banner; banner = str)

function print_banner()
    println(banner)
end

function print_command_help()
    print("""
        usage: $name command options...

        commands:
        """)
    for command in commands
        println(command)
    end
    println("")
end

function print_option_help(command)
    command_options = options[command]
    println("$command options:")
    for option in command_options
        println(option)
    end
    println("")
end

function parse_option(opt::Option,args)
    if length(args) > opt.max
        error("Too many arguments passed to flag $(opt.flag) ($(opt.max) maximum).")
    end
    if length(args) < opt.min
        error("Too few arguments passed to flag $(opt.flag) ($(opt.min) minimum).")
    end
    if opt.max == 0
        return nothing
    end
    if opt.max > 1
        return opt.T[parse_option(opt.T,arg) for arg in args]
    end
    parse_option(opt.T,args[1])
end

function parse_option(T::Type,arg)
    try
        return parse_option_helper(T,arg)
    catch
        error("Unable to convert $arg to type $T.")
    end
end

parse_option_helper(::Type{UTF8String},arg) = utf8(arg)
parse_option_helper(::Type{Int},arg) = parse(Int,arg)
parse_option_helper(::Type{Float64},arg) = parse(Float64,arg)

looks_like_flag(str) = startswith(str,"--")

function parse_args(args)
    command_names = [command.name for command in commands]
    if length(args) < 1 || !(args[1] in command_names)
        print_command_help()
        error("Please provide one of the listed commands.")
    end

    command = args[1]
    option_flags = [option.flag for option in options[command]]

    if length(args) < 2
        print_option_help(command)
        error("Please select from the list options.")
    end

    args_dict = Dict{UTF8String,Any}()
    idx = 2
    while idx <= length(args)
        if !(args[idx] in option_flags)
            print_option_help(command)
            error("\"$(args[idx])\" is not a recognized flag.")
        end

        option_flag = args[idx]
        option = options[command][findfirst(option_flags,option_flag)]

        next_idx = findnext(looks_like_flag,args,idx+1)
        if next_idx == 0 # This is the last option
            args_dict[option_flag] = parse_option(option,args[idx+1:end])
            break
        else
            args_dict[option_flag] = parse_option(option,args[idx+1:next_idx-1])
            idx = next_idx
        end
    end

    # Verify that the required options are provided
    for option in options[command]
        if option.required && !(option.flag in keys(args_dict))
            error("Required flag $(option.flag) is missing.")
        end
    end

    # Verify that there are no conflicts
    for option_flag in keys(args_dict)
        option = options[command][findfirst(option_flags,option_flag)]
        if !isempty(intersect(option.conflicts,keys(args_dict)))
            error("There is a conflict with option $(option_flag).")
        end
    end

    command,args_dict
end

end

