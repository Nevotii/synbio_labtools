# -*- coding: utf-8 -*-

from collections import OrderedDict
import streamlit as st
from streamlit.logger import get_logger
import tools


LOGGER = get_logger(__name__)
# Dictionary of
# demo_name -> (demo_function, demo_description)

TOOLS = OrderedDict(
    [("—", (tools.intro, None)),
        ("Gibson calculator",
            (tools.gibson,
                """A small gibson calculator for quickly organising your reaction""",
            ),
        ),
        ("C. acnes genome pattern matching",
            (tools.pattern,
                """To quickly find arbitrary sequences in C. acnes genome.""",
            ),
        )     ,
        ("Base editor finder",
            (tools.baseedit,
                """To quickly find arbitrary sequences in C. acnes genome.""",
            ),
        )
    ]
)

def run():
    tool_name = st.sidebar.selectbox("Choose a tool", list(TOOLS.keys()), 0)
    tool = TOOLS[tool_name][0]

    if tool_name =='—':
        st.write("# Welcome to Synbio Tools 👋")

    tool()


    
    for i in range(10):
        st.empty()
        




    st.markdown("***")
    st.info("This page was created by [Guillermo Nevot](https://gnevot.xyz). If you like the app, please [reach out!](https://twitter.com/nevotii?lang=en)")


if __name__ == "__main__":
    run()
