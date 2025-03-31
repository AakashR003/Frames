class Config:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(Config, cls).__new__(cls)
            cls._instance.FEDivision = 20  # Default value
        return cls._instance

    def get_FEDivision(self):
        return self.FEDivision  # Always return the latest value

    def set_FEDivision(self, value):  # ✅ New setter function
        self.FEDivision = value  # Update the value dynamically

config = Config()  # Singleton instance
