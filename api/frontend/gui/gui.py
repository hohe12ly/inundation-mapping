import os
from gevent import monkey
from flask import Flask, render_template, request

monkey.patch_all()

NODE_GUI_CONNECTOR_URL = os.environ.get('NODE_GUI_CONNECTOR_URL')

app = Flask(__name__)

@app.route('/')
def main():
    return render_template('index.html', node_gui_connector_url=NODE_GUI_CONNECTOR_URL)

if __name__ == '__main__':
    app.run("0.0.0.0", port=5000)
