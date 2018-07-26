
set -x
set -e
ulimit -s unlimited
ulimit -c 0

function error_exit
{
  if [ $1 -ne 0 ]; then
    echo "Error with exit code ${1}"
    if [ -e FrameworkJobReport.xml ]; then
      cat << EOF > FrameworkJobReport.xml.tmp
      <FrameworkJobReport>
      <FrameworkError ExitStatus="${1}" Type="" >
      Error with exit code ${1}
      </FrameworkError>
EOF
      tail -n+2 FrameworkJobReport.xml >> FrameworkJobReport.xml.tmp
      mv FrameworkJobReport.xml.tmp FrameworkJobReport.xml
    else
      cat << EOF > FrameworkJobReport.xml
      <FrameworkJobReport>
      <FrameworkError ExitStatus="${1}" Type="" >
      Error with exit code ${1}
      </FrameworkError>
      </FrameworkJobReport>
EOF
    fi
    exit 0
  fi
}

trap 'error_exit $?' ERR

while IFS='' read -r line || [[ -n "$line" ]]; do
  filesarray+=( "$line" )
done < VLQAnalysis_Ntuples.txt

if [ $1 -eq 1 ]; then
  python testVLQAnalyzer.py -f ${filesarray[0]}
fi
if [ $1 -eq 2 ]; then
  python testVLQAnalyzer.py -f ${filesarray[1]}
fi
if [ $1 -eq 3 ]; then
  python testVLQAnalyzer.py -f ${filesarray[2]}
fi
if [ $1 -eq 4 ]; then
  python testVLQAnalyzer.py -f ${filesarray[3]}
fi
if [ $1 -eq 5 ]; then
  python testVLQAnalyzer.py -f ${filesarray[4]}
fi
if [ $1 -eq 6 ]; then
  python testVLQAnalyzer.py -f ${filesarray[5]}
fi
if [ $1 -eq 7 ]; then
  python testVLQAnalyzer.py -f ${filesarray[6]}
fi
if [ $1 -eq 8 ]; then
  python testVLQAnalyzer.py -f ${filesarray[7]}
fi
if [ $1 -eq 9 ]; then
  python testVLQAnalyzer.py -f ${filesarray[8]}
fi
if [ $1 -eq 10 ]; then
  python testVLQAnalyzer.py -f ${filesarray[9]}
fi
if [ $1 -eq 11 ]; then
  python testVLQAnalyzer.py -f ${filesarray[10]}
fi
if [ $1 -eq 12 ]; then
  python testVLQAnalyzer.py -f ${filesarray[11]}
fi
if [ $1 -eq 13 ]; then
  python testVLQAnalyzer.py -f ${filesarray[12]}
fi
if [ $1 -eq 14 ]; then
  python testVLQAnalyzer.py -f ${filesarray[13]}
fi
if [ $1 -eq 15 ]; then
  python testVLQAnalyzer.py -f ${filesarray[14]}
fi
if [ $1 -eq 16 ]; then
  python testVLQAnalyzer.py -f ${filesarray[15]}
fi
if [ $1 -eq 17 ]; then
  python testVLQAnalyzer.py -f ${filesarray[16]}
fi
if [ $1 -eq 18 ]; then
  python testVLQAnalyzer.py -f ${filesarray[17]}
fi
if [ $1 -eq 19 ]; then
  python testVLQAnalyzer.py -f ${filesarray[18]}
fi
if [ $1 -eq 20 ]; then
  python testVLQAnalyzer.py -f ${filesarray[19]}
fi
if [ $1 -eq 21 ]; then
  python testVLQAnalyzer.py -f ${filesarray[20]}
fi
if [ $1 -eq 22 ]; then
  python testVLQAnalyzer.py -f ${filesarray[21]}
fi
if [ $1 -eq 23 ]; then
  python testVLQAnalyzer.py -f ${filesarray[22]}
fi
if [ $1 -eq 24 ]; then
  python testVLQAnalyzer.py -f ${filesarray[23]}
fi
if [ $1 -eq 25 ]; then
  python testVLQAnalyzer.py -f ${filesarray[24]}
fi
if [ $1 -eq 26 ]; then
  python testVLQAnalyzer.py -f ${filesarray[25]}
fi
if [ $1 -eq 27 ]; then
  python testVLQAnalyzer.py -f ${filesarray[26]}
fi
if [ $1 -eq 28 ]; then
  python testVLQAnalyzer.py -f ${filesarray[27]}
fi
if [ $1 -eq 29 ]; then
  python testVLQAnalyzer.py -f ${filesarray[28]}
fi
tar -cf VLQAnalysis_output.tar *out.root
rm *out.root
