import os

from flask import Flask, render_template, url_for, request, send_file, flash

from motif_search import match_list_from_uniprot

# Create the application object:
app = Flask(__name__)
# Generate secret key (necessary for flash messages since they use a cookie):
app.config['SECRET_KEY'] = os.urandom(16)
# Use decorators to link the function to a url:
@app.route("/", methods=['GET', 'POST'])
def home():
    """Home page function, takes in form input and renders home.html"""
    # Post method is called if submit button is pressed.
    if request.method == 'POST':
        # Get the contents of the "uniprot" text field:
        ac_code = request.form["uniprot"]
        # Get processed results, generate csv file for download:
        results = match_list_from_uniprot(ac_code)
        # If ac_code is invalid, results will be None:
        if results is None:
            flash("No protein record found, please enter valid accession number.", "danger")
        # Set uniprot url of protein:
        prot_url = "https://uniprot.org/uniprot/" + ac_code
    else:
        # Set all passed variables to None to avoid triggering the if blocks in html:
        results, ac_code, prot_url = None, None, None
    # Render home.html, passing in the results, accession code, and uniprot url:
    return render_template('home.html', results=results, ac_code=ac_code, prot_url=prot_url)

# Route for downloading csv file:
@app.route('/return-file/')
def return_file():
    """This page exists just to download files, it will open, download file, and
    immediately autoclose when the download button is pressed."""
    try:
        # This gets results.csv from the home directory and downloads it:
        return send_file('results.csv', as_attachment=True)
    except Exception as e:
        return str(e)

# Create about page:
@app.route("/about")
def about():
    return render_template('about.html', title='About')

# Start the server with the 'run()' method:
if __name__ == '__main__':
    app.run()
