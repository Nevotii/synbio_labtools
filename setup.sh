mkdir -p ~/.streamlit/
rm ~/.streamlit/config.toml
echo "\
[general]\n\
email = \"guillermo.nevot@upf.edu\"\n\
" > ~/.streamlit/credentials.toml
echo "\
[server]\n\
headless = true\n\
enableCORS=false\n\
port = $PORT\n\
" > ~/.streamlit/config.toml
